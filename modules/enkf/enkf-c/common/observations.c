/******************************************************************************
 *
 * File:        observations.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdint.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "observations.h"

#define NOBSTYPES_INC 10
#define KD_INC 50000
#define HT_SIZE 500

typedef struct {
    int obstypeid;
    char fname[MAXSTRLEN];
    int fid;
    int batch;
} badbatch;

typedef union {
    int key_int[2];
    short int key_short[4];
} keydata;

/**
 */
void obs_addtype(observations* obs, obstype* src)
{
    obstype* ot;
    int i;

    if (obs->nobstypes % NOBSTYPES_INC == 0)
        obs->obstypes = realloc(obs->obstypes, (obs->nobstypes + NOBSTYPES_INC) * sizeof(obstype));
    ot = &obs->obstypes[obs->nobstypes];
    ot->id = obs->nobstypes;
    ot->name = strdup(src->name);
    ot->issurface = src->issurface;
    ot->nvar = src->nvar;
    ot->varnames = malloc(src->nvar * sizeof(char*));
    for (i = 0; i < src->nvar; ++i)
        ot->varnames[i] = strdup(src->varnames[i]);
    ot->alias = strdup(src->alias);
    ot->offset_fname = (src->offset_fname != NULL) ? strdup(src->offset_fname) : NULL;
    ot->offset_varname = (src->offset_varname != NULL) ? strdup(src->offset_varname) : NULL;
    ot->mld_varname = (src->mld_varname != NULL) ? strdup(src->mld_varname) : NULL;
    ot->mld_threshold = src->mld_threshold;
    ot->hfunction = strdup(src->hfunction);
    ot->allowed_min = src->allowed_min;
    ot->allowed_max = src->allowed_max;
    ot->isasync = src->isasync;
    ot->async_tstep = src->async_tstep;
    ot->nlocrad = src->nlocrad;
    ot->locrad = malloc(sizeof(double) * ot->nlocrad);
    ot->locweight = malloc(sizeof(double) * ot->nlocrad);
    for (i = 0; i < ot->nlocrad; ++i) {
        ot->locrad[i] = src->locrad[i];
        ot->locweight[i] = src->locweight[i];
    }
    ot->rfactor = src->rfactor;
    ot->nlobsmax = src->nlobsmax;
    ot->estdmin = src->estdmin;
    ot->nsubgrid = src->nsubgrid;
    ot->nmodified = 0;
    ot->date_min = src->date_min;
    ot->date_max = src->date_max;
    ot->xmin = src->xmin;
    ot->xmax = src->xmax;
    ot->ymin = src->ymin;
    ot->ymax = src->ymax;
    ot->zmin = src->zmin;
    ot->zmax = src->zmax;

    ot->ndomains = src->ndomains;
    if (ot->ndomains > 0)
        ot->domainnames = malloc(ot->ndomains * sizeof(char*));
    for (i = 0; i < ot->ndomains; ++i)
        ot->domainnames[i] = strdup(src->domainnames[i]);

    /*
     * (these fields are set by obs_calcstats())
     */
    ot->nobs = -1;
    ot->ngood = -1;
    ot->noutside_grid = -1;
    ot->noutside_obsdomain = -1;
    ot->noutside_obswindow = -1;
    ot->nland = -1;
    ot->nshallow = -1;
    ot->nbadbatch = -1;
    ot->nrange = -1;
    ot->nthinned = -1;

    obs->nobstypes++;
}

/**
 */
observations* obs_create(void)
{
    observations* obs = malloc(sizeof(observations));

    obs->products = NULL;
    obs->instruments = NULL;
    obs->datafiles = NULL;
    obs->nobstypes = 0;
    obs->obstypes = NULL;
#if defined(ENKF_CALC)
    obs->loctrees = NULL;
#endif
    obs->obsids = NULL;
    obs->da_date = NAN;
    obs->datestr = NULL;
    obs->allobs = 0;
    obs->nallocated = 0;
    obs->nobs = 0;
    obs->data = NULL;
    obs->compacted = 0;
    obs->hasstats = 0;
    obs->ngood = 0;
    obs->noutside_grid = 0;
    obs->noutside_obsdomain = 0;
    obs->noutside_obswindow = 0;
    obs->nland = 0;
    obs->nshallow = 0;
    obs->nthinned = 0;
    obs->nbadbatch = 0;
    obs->nrange = 0;
    obs->nmodified = 0;
    obs->badbatches = NULL;
    obs->ncformat = NETCDF_FORMAT;
    obs->nccompression = 0;
#if defined(ENKF_PREP)
    obs->model = NULL;
#endif

    return obs;
}

/**
 */
observations* obs_create_fromprm(enkfprm* prm)
{
    observations* obs = obs_create();

    obs->products = st_create("products");
    obs->instruments = st_create("instruments");
    obs->datafiles = st_create("datafiles");

    enkf_printf("  reading observation type specs from \"%s\":\n", prm->obstypeprm);
    obstypes_read(prm, prm->obstypeprm, &obs->nobstypes, &obs->obstypes);

#if defined(ENKF_PREP)
    obs->da_date = date_str2dbl(prm->date);
    obs->datestr = strdup(prm->date);

    if (file_exists(FNAME_BADBATCHES)) {
        FILE* f = NULL;
        char buf[MAXSTRLEN];
        int line;

        enkf_printf("  reading bad batches:\n");

        obs->badbatches = ht_create_i1s2(HT_SIZE);
        f = enkf_fopen(FNAME_BADBATCHES, "r");
        line = 0;
        while (fgets(buf, MAXSTRLEN, f) != NULL) {
            char obstype[MAXSTRLEN];
            badbatch* bb;
            keydata key;

            line++;
            if (buf[0] == '#')
                continue;
            bb = malloc(sizeof(badbatch));
            if (sscanf(buf, "%s %s %d %d", obstype, bb->fname, &bb->fid, &bb->batch) != 4)
                 enkf_quit("%s, l.%d: wrong bad batch specification (expected \"%s %s %d %d\"\n", FNAME_BADBATCHES, line);

            bb->obstypeid = obstype_getid(obs->nobstypes, obs->obstypes, obstype, 1);

            key.key_int[0] = bb->batch;
            key.key_short[2] = bb->obstypeid;
            key.key_short[3] = bb->fid;
            ht_insert(obs->badbatches, &key, bb);
            enkf_printf("    %s %s %d %d\n", obstype, bb->fname, bb->fid, bb->batch);
        }
        fclose(f);
    }

    {
        int otid;

        for (otid = 0; otid <= obs->nobstypes; ++otid) {
            obstype* ot = &obs->obstypes[otid];

            if (isnan(ot->windowmin))
                ot->windowmin = prm->windowmin;
            if (isnan(ot->windowmax))
                ot->windowmax = prm->windowmax;
        }
    }
#endif

    obs->ncformat = prm->ncformat;
    obs->nccompression = prm->nccompression;

    return obs;
}

#if defined(ENKF_PREP)
/**
 */
#define BYTE_PER_SHORT 2
#define BYTE_PER_INT 4
void obs_markbadbatches(observations* obs)
{
    int i;

    assert(BYTE_PER_SHORT == sizeof(short));
    assert(BYTE_PER_INT == sizeof(int));

    if (obs->badbatches == NULL || ht_getnentries(obs->badbatches) == 0)
        return;

    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];
        keydata key;
        badbatch* bb;

        key.key_int[0] = o->batch;
        key.key_short[2] = o->type;
        key.key_short[3] = o->fid;
        bb = ht_find(obs->badbatches, &key);

        if (bb != NULL && bb->obstypeid == o->type) {
            if (strcmp(bb->fname, st_findstringbyindex(obs->datafiles, o->fid)) != 0)
                enkf_quit("bad batch processing: file name for fid = %d in \"%s\" does not match the data file name. Check that \"%s\" is the right file.", o->fid, FNAME_BADBATCHES, FNAME_BADBATCHES);
            o->status = STATUS_BADBATCH;
        }
    }
}
#endif

/**
 */
observations* obs_create_fromdata(observations* parentobs, int nobs, observation data[])
{
    observations* obs = obs_create();
    int i;

    obs->products = st_copy(parentobs->products);
    obs->instruments = st_copy(parentobs->instruments);
    obs->datafiles = st_copy(parentobs->datafiles);

    for (i = 0; i < parentobs->nobstypes; ++i)
        obs_addtype(obs, &parentobs->obstypes[i]);

    obs->da_date = parentobs->da_date;
    obs->datestr = strdup(parentobs->datestr);

    obs->nobs = nobs;
    obs->data = data;

    obs->ncformat = parentobs->ncformat;
    obs->nccompression = parentobs->nccompression;

    return obs;
}

/**
 */
void obs_destroy(observations* obs)
{
    st_destroy(obs->products);
    st_destroy(obs->instruments);
    st_destroy(obs->datafiles);
    obstypes_destroy(obs->nobstypes, obs->obstypes);
#if defined(ENKF_CALC)
    if (obs->loctrees != NULL) {
        int i;

        for (i = 0; i < obs->nobstypes; ++i)
            if (obs->loctrees[i] != NULL)
                kd_destroy(obs->loctrees[i]);
        free(obs->loctrees);
    }
#endif
    if (obs->obsids != NULL) {
        int i;

        for (i = 0; i < obs->nobstypes; ++i)
            if (obs->obsids[i] != NULL)
                free(obs->obsids[i]);
        free(obs->obsids);
    }
    if (obs->nobs > 0)
        free(obs->data);
    if (obs->datestr != NULL)
        free(obs->datestr);
    if (obs->badbatches != NULL) {
        ht_process(obs->badbatches, free);
        ht_destroy(obs->badbatches);
    }
    free(obs);
}

/**
 */
void obs_checkalloc(observations* obs)
{
    if (obs->nobs == obs->nallocated) {
        obs->data = realloc(obs->data, (obs->nobs + NOBS_INC) * sizeof(observation));
        if (obs->data == NULL)
            enkf_quit("obs_checkalloc(): not enough memory");
        obs->nallocated += NOBS_INC;
        memset(&obs->data[obs->nobs], 0, NOBS_INC * sizeof(observation));
    }
}

/**
 */
static int comp_obsstatus(const void* p1, const void* p2)
{
    observation* m1 = (observation*) p1;
    observation* m2 = (observation*) p2;

    if (m1->status > m2->status)
        return 1;
    if (m1->status < m2->status)
        return -1;
    return 0;
}

/**
 */
void obs_compact(observations* obs)
{
    int i;

    if (obs->compacted)
        return;

    enkf_printf("    compacting obs:");
    enkf_flush();
    assert(STATUS_OK == 0);
    qsort(obs->data, obs->nobs, sizeof(observation), comp_obsstatus);
    for (i = 0; i < obs->nobs; ++i) {
        obs->data[i].id_orig = obs->data[i].id;
        obs->data[i].id = i;
    }
    enkf_printf("\n");
    enkf_flush();
    obs->compacted = 1;
}

/**
 */
static int comp_obsid(const void* p1, const void* p2)
{
    observation* m1 = (observation*) p1;
    observation* m2 = (observation*) p2;

    if (m1->id > m2->id)
        return 1;
    if (m1->id < m2->id)
        return -1;
    return 0;
}

/** Sort observations by id.
 */
void obs_inorder(observations* obs)
{
    qsort(obs->data, obs->nobs, sizeof(observation), comp_obsid);
}

/**
 */
void obs_calcstats(observations* obs)
{
    int i;

    if (obs->hasstats)
        return;

    obs->ngood = 0;
    obs->noutside_grid = 0;
    obs->noutside_obsdomain = 0;
    obs->noutside_obswindow = 0;
    obs->nland = 0;
    obs->nshallow = 0;
    obs->nrange = 0;
    for (i = 0; i < obs->nobstypes; ++i) {
        obstype* ot = &obs->obstypes[i];

        ot->nobs = 0;
        ot->ngood = 0;
        ot->noutside_grid = 0;
        ot->noutside_obsdomain = 0;
        ot->noutside_obswindow = 0;
        ot->nland = 0;
        ot->nshallow = 0;
        ot->nbadbatch = 0;
        ot->nrange = 0;
        ot->date_min = DBL_MAX;
        ot->date_max = -DBL_MAX;
    }

    for (i = 0; i < obs->nobs; ++i) {
        observation* m = &obs->data[i];
        obstype* ot = &obs->obstypes[m->type];

        ot->nobs++;
        if (m->status == STATUS_OK) {
            obs->ngood++;
            ot->ngood++;
        } else if (m->status == STATUS_OUTSIDEGRID) {
            obs->noutside_grid++;
            ot->noutside_grid++;
        } else if (m->status == STATUS_OUTSIDEOBSDOMAIN) {
            obs->noutside_obsdomain++;
            ot->noutside_obsdomain++;
        } else if (m->status == STATUS_OUTSIDEOBSWINDOW) {
            obs->noutside_obswindow++;
            ot->noutside_obswindow++;
        } else if (m->status == STATUS_LAND) {
            obs->nland++;
            ot->nland++;
        } else if (m->status == STATUS_SHALLOW) {
            obs->nshallow++;
            ot->nshallow++;
        } else if (m->status == STATUS_BADBATCH) {
            obs->nbadbatch++;
            ot->nbadbatch++;
        } else if (m->status == STATUS_RANGE) {
            obs->nrange++;
            ot->nrange++;
        } else if (m->status == STATUS_THINNED) {
            obs->nthinned++;
            ot->nthinned++;
        }

        if (m->date < ot->date_min)
            ot->date_min = m->date;
        if (m->date > ot->date_max)
            ot->date_max = m->date;
    }
    obs->hasstats = 1;
}

/** Reads observations from "observations.nc".
 */
void obs_read(observations* obs, char fname[])
{
    int ncid;
    double da_julday = NAN;
    int dimid_nobs[1];
    size_t nobs;
    int varid_type, varid_product, varid_instrument, varid_id, varid_idorig, varid_fid, varid_batch, varid_value, varid_std, varid_lon, varid_lat, varid_depth, varid_mdepth, varid_fi, varid_fj, varid_fk, varid_date, varid_status, varid_aux;
    int* id;
    int* id_orig;
    short int* type;
    short int* product;
    short int* instrument;
    short int* fid;
    int* batch;
    double* value;
    double* std;
    double* lon;
    double* lat;
    double* depth;
    double* model_depth;
    double* fi;
    double* fj;
    double* fk;
    double* date;
    int* status;
    int* aux;
    int natts;
    int i;

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_get_att_double(ncid, NC_GLOBAL, "DA_JULDAY", &da_julday);
    if (!enkf_noobsdatecheck && (isnan(da_julday) || fabs(obs->da_date - da_julday) > 1e-6))
        enkf_quit("\"observations.nc\" from a different cycle");

    ncw_inq_dimid(ncid, "nobs", dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs[0], &nobs);

    obs->nobs = nobs;
    enkf_printf("    %zu observations\n", nobs);
    if (nobs == 0) {
        obs->data = NULL;
        ncw_close(ncid);
        goto finish;
    }

    enkf_printf("    allocating %zu bytes for array of observations\n", nobs * sizeof(observation));
    obs->data = malloc(nobs * sizeof(observation));
    assert(obs->data != NULL);
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    ncw_inq_varid(ncid, "type", &varid_type);
    ncw_inq_varid(ncid, "product", &varid_product);
    ncw_inq_varid(ncid, "instrument", &varid_instrument);
    ncw_inq_varid(ncid, "id", &varid_id);
    ncw_inq_varid(ncid, "id_orig", &varid_idorig);
    ncw_inq_varid(ncid, "fid", &varid_fid);
    ncw_inq_varid(ncid, "batch", &varid_batch);
    ncw_inq_varid(ncid, "value", &varid_value);
    ncw_inq_varid(ncid, "std", &varid_std);
    ncw_inq_varid(ncid, "lon", &varid_lon);
    ncw_inq_varid(ncid, "lat", &varid_lat);
    ncw_inq_varid(ncid, "depth", &varid_depth);
    ncw_inq_varid(ncid, "model_depth", &varid_mdepth);
    ncw_inq_varid(ncid, "fi", &varid_fi);
    ncw_inq_varid(ncid, "fj", &varid_fj);
    ncw_inq_varid(ncid, "fk", &varid_fk);
    ncw_inq_varid(ncid, "date", &varid_date);
    ncw_inq_varid(ncid, "status", &varid_status);
    ncw_inq_varid(ncid, "aux", &varid_aux);

    /*
     * date
     */
    {
        size_t len = 0;

        ncw_inq_attlen(ncid, varid_date, "units", &len);
        obs->datestr = malloc(len + 1);
        ncw_get_att_text(ncid, varid_date, "units", obs->datestr);
        obs->datestr[len] = 0;
    }

    /*
     * type (basically, just a check)
     */
    ncw_inq_varnatts(ncid, varid_type, &natts);
    for (i = 0; i < natts; ++i) {
        char attname[NC_MAX_NAME];
        int typeid;

        ncw_inq_attname(ncid, varid_type, i, attname);
        typeid = obstype_getid(obs->nobstypes, obs->obstypes, attname, 0);
        if (typeid >= 0) {
            int typeid_read;

            ncw_check_attlen(ncid, varid_type, attname, 1);
            ncw_get_att_int(ncid, varid_type, attname, &typeid_read);
            assert(typeid == typeid_read);
        }
    }

    /*
     * product 
     */
    ncw_inq_varnatts(ncid, varid_product, &natts);
    for (i = 0; i < natts; ++i) {
        char name[NC_MAX_NAME];
        nc_type type;
        size_t len;

        ncw_inq_attname(ncid, varid_product, i, name);
        ncw_inq_att(ncid, varid_product, name, &type, &len);
        if (type == NC_INT && len == 1) {
            int id;

            ncw_get_att_int(ncid, varid_product, name, &id);
            st_add(obs->products, name, id);
        }
    }

    /*
     * instrument 
     */
    ncw_inq_varnatts(ncid, varid_instrument, &natts);
    for (i = 0; i < natts; ++i) {
        char name[NC_MAX_NAME];
        nc_type type;
        size_t len;

        ncw_inq_attname(ncid, varid_instrument, i, name);
        ncw_inq_att(ncid, varid_instrument, name, &type, &len);
        if (type == NC_INT && len == 1) {
            int id;

            ncw_get_att_int(ncid, varid_instrument, name, &id);
            st_add(obs->instruments, name, id);
        }
    }

    /*
     * datafiles
     */
    ncw_inq_varnatts(ncid, varid_fid, &natts);
    for (i = 0; i < natts; ++i) {
        char name[NC_MAX_NAME];
        char attstr[MAXSTRLEN];
        size_t len;
        int id;

        ncw_inq_attname(ncid, varid_fid, i, name);
        if (!str2int(name, &id))
            continue;
        ncw_inq_attlen(ncid, varid_fid, name, &len);
        assert(len < MAXSTRLEN);
        ncw_get_att_text(ncid, varid_fid, name, attstr);
        attstr[len] = 0;
        st_add_ifabsent(obs->datafiles, attstr, id);
    }

    id = malloc(nobs * sizeof(int));
    id_orig = malloc(nobs * sizeof(int));
    type = malloc(nobs * sizeof(short int));
    product = malloc(nobs * sizeof(short int));
    instrument = malloc(nobs * sizeof(short int));
    fid = malloc(nobs * sizeof(short int));
    batch = malloc(nobs * sizeof(int));
    value = malloc(nobs * sizeof(double));
    std = malloc(nobs * sizeof(double));
    lon = malloc(nobs * sizeof(double));
    lat = malloc(nobs * sizeof(double));
    depth = malloc(nobs * sizeof(double));
    model_depth = malloc(nobs * sizeof(double));
    fi = malloc(nobs * sizeof(double));
    fj = malloc(nobs * sizeof(double));
    fk = malloc(nobs * sizeof(double));
    date = malloc(nobs * sizeof(double));
    status = malloc(nobs * sizeof(int));
    aux = malloc(nobs * sizeof(int));

    ncw_get_var_int(ncid, varid_id, id);
    ncw_get_var_int(ncid, varid_idorig, id_orig);
    ncw_get_var_short(ncid, varid_type, type);
    ncw_get_var_short(ncid, varid_product, product);
    ncw_get_var_short(ncid, varid_instrument, instrument);
    ncw_get_var_short(ncid, varid_fid, fid);
    ncw_get_var_int(ncid, varid_batch, batch);
    ncw_get_var_double(ncid, varid_value, value);
    ncw_get_var_double(ncid, varid_std, std);
    ncw_get_var_double(ncid, varid_lon, lon);
    ncw_get_var_double(ncid, varid_lat, lat);
    ncw_get_var_double(ncid, varid_depth, depth);
    ncw_get_var_double(ncid, varid_mdepth, model_depth);
    ncw_get_var_double(ncid, varid_fi, fi);
    ncw_get_var_double(ncid, varid_fj, fj);
    ncw_get_var_double(ncid, varid_fk, fk);
    ncw_get_var_double(ncid, varid_date, date);
    ncw_get_var_int(ncid, varid_status, status);
    ncw_get_var_int(ncid, varid_aux, aux);

    ncw_close(ncid);

    for (i = 0; i < (int) nobs; ++i) {
        observation* m = &obs->data[i];

        m->type = type[i];
        m->product = product[i];
        m->instrument = instrument[i];
        m->id = id[i];
        m->id_orig = id_orig[i];
        m->fid = fid[i];
        m->batch = batch[i];
        m->value = value[i];
        m->std = std[i];
        m->lon = lon[i];
        m->lat = lat[i];
        m->depth = depth[i];
        m->model_depth = model_depth[i];
        m->fi = fi[i];
        m->fj = fj[i];
        m->fk = fk[i];
        m->date = date[i];
        m->status = status[i];
        m->aux = aux[i];
    }

    free(type);
    free(product);
    free(instrument);
    free(id);
    free(id_orig);
    free(fid);
    free(batch);
    free(value);
    free(std);
    free(lon);
    free(lat);
    free(depth);
    free(model_depth);
    free(fi);
    free(fj);
    free(fk);
    free(date);
    free(status);
    free(aux);
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

  finish:

    obs_calcstats(obs);
}

/**
 */
void obs_write(observations* obs, char fname[])
{
    int nobs = obs->nobs;
    char tunits[MAXSTRLEN];

    int ncid;
    int dimid_nobs[1];
    int varid_type, varid_product, varid_instrument, varid_id, varid_idorig, varid_fid, varid_batch, varid_value, varid_std, varid_lon, varid_lat, varid_depth, varid_mdepth, varid_fi, varid_fj, varid_fk, varid_date, varid_status, varid_aux;

    int* id;
    int* id_orig;
    short int* type;
    short int* product;
    short int* instrument;
    short int* fid;
    int* batch;
    double* value;
    double* std;
    double* lon;
    double* lat;
    double* depth;
    double* model_depth;
    double* fi;
    double* fj;
    double* fk;
    double* date;
    int* status;
    int* aux;

    int i, ii;

    if (rank != 0)
        return;

    if (file_exists(fname))
        enkf_quit("file \"%s\" already exists", fname);
    ncw_create(fname, NC_NOCLOBBER | obs->ncformat, &ncid);

    ncw_put_att_double(ncid, NC_GLOBAL, "DA_JULDAY", 1, &obs->da_date);

    ncw_def_dim(ncid, "nobs", nobs, dimid_nobs);
    ncw_def_var(ncid, "id", NC_INT, 1, dimid_nobs, &varid_id);
    ncw_put_att_text(ncid, varid_id, "long_name", "observation ID");
    ncw_def_var(ncid, "id_orig", NC_INT, 1, dimid_nobs, &varid_idorig);
    ncw_put_att_text(ncid, varid_idorig, "long_name", "original observation ID");
    ncw_put_att_text(ncid, varid_idorig, "description", "for primary observations - the serial number of the primary observation during the reading of data files; for superobs - the original ID of the very first observation collated into this observation");
    ncw_def_var(ncid, "type", NC_SHORT, 1, dimid_nobs, &varid_type);
    ncw_put_att_text(ncid, varid_type, "long_name", "observation type ID");
    ncw_def_var(ncid, "product", NC_SHORT, 1, dimid_nobs, &varid_product);
    ncw_put_att_text(ncid, varid_product, "long_name", "observation product ID");
    ncw_def_var(ncid, "instrument", NC_SHORT, 1, dimid_nobs, &varid_instrument);
    ncw_put_att_text(ncid, varid_instrument, "long_name", "observation instrument ID");
    ncw_def_var(ncid, "fid", NC_SHORT, 1, dimid_nobs, &varid_fid);
    ncw_put_att_text(ncid, varid_fid, "long_name", "observation data file ID");
    ncw_def_var(ncid, "batch", NC_INT, 1, dimid_nobs, &varid_batch);
    ncw_put_att_text(ncid, varid_batch, "long_name", "observation batch ID");
    ncw_def_var(ncid, "value", NC_FLOAT, 1, dimid_nobs, &varid_value);
    ncw_put_att_text(ncid, varid_value, "long_name", "observation value");
    ncw_def_var(ncid, "std", NC_FLOAT, 1, dimid_nobs, &varid_std);
    ncw_put_att_text(ncid, varid_std, "long_name", "standard deviation of observation error used in DA");
    ncw_def_var(ncid, "lon", NC_FLOAT, 1, dimid_nobs, &varid_lon);
    ncw_put_att_text(ncid, varid_lon, "long_name", "observation longitude");
    ncw_def_var(ncid, "lat", NC_FLOAT, 1, dimid_nobs, &varid_lat);
    ncw_put_att_text(ncid, varid_lat, "long_name", "observation latitude");
    ncw_def_var(ncid, "depth", NC_FLOAT, 1, dimid_nobs, &varid_depth);
    ncw_put_att_text(ncid, varid_depth, "long_name", "observation depth/height");
    ncw_def_var(ncid, "model_depth", NC_FLOAT, 1, dimid_nobs, &varid_mdepth);
    ncw_put_att_text(ncid, varid_mdepth, "long_name", "model bottom depth at the observation location");
    ncw_def_var(ncid, "fi", NC_FLOAT, 1, dimid_nobs, &varid_fi);
    ncw_put_att_text(ncid, varid_fi, "long_name", "fractional grid index i of the observation");
    ncw_def_var(ncid, "fj", NC_FLOAT, 1, dimid_nobs, &varid_fj);
    ncw_put_att_text(ncid, varid_fj, "long_name", "fractional grid index j of the observation");
    ncw_def_var(ncid, "fk", NC_FLOAT, 1, dimid_nobs, &varid_fk);
    ncw_put_att_text(ncid, varid_fk, "long_name", "fractional grid index k of the observation");
    ncw_def_var(ncid, "date", NC_FLOAT, 1, dimid_nobs, &varid_date);
    ncw_put_att_text(ncid, varid_date, "long_name", "observation time");
    ncw_def_var(ncid, "status", NC_BYTE, 1, dimid_nobs, &varid_status);
    ncw_put_att_text(ncid, varid_status, "long_name", "observation status");
    i = STATUS_OK;
    ncw_put_att_int(ncid, varid_status, "STATUS_OK", 1, &i);
    i = STATUS_OUTSIDEGRID;
    ncw_put_att_int(ncid, varid_status, "STATUS_OUTSIDEGRID", 1, &i);
    i = STATUS_LAND;
    ncw_put_att_int(ncid, varid_status, "STATUS_LAND", 1, &i);
    i = STATUS_SHALLOW;
    ncw_put_att_int(ncid, varid_status, "STATUS_SHALLOW", 1, &i);
    i = STATUS_RANGE;
    ncw_put_att_int(ncid, varid_status, "STATUS_RANGE", 1, &i);
    i = STATUS_BADBATCH;
    ncw_put_att_int(ncid, varid_status, "STATUS_BADBATCH", 1, &i);
    i = STATUS_OUTSIDEOBSDOMAIN;
    ncw_put_att_int(ncid, varid_status, "STATUS_OUTSIDEOBSDOMAIN", 1, &i);
    i = STATUS_OUTSIDEOBSWINDOW;
    ncw_put_att_int(ncid, varid_status, "STATUS_OUTSIDEOBSWINDOW", 1, &i);
    i = STATUS_THINNED;
    ncw_put_att_int(ncid, varid_status, "STATUS_THINNED", 1, &i);
    ncw_def_var(ncid, "aux", NC_INT, 1, dimid_nobs, &varid_aux);
    ncw_put_att_text(ncid, varid_aux, "long_name", "auxiliary information");
    ncw_put_att_text(ncid, varid_aux, "description", "for primary observations - the ID of the superobservation it is collated into; for superobservations - the number of primary observations collated");
    snprintf(tunits, MAXSTRLEN, "days from %s", obs->datestr);
    ncw_put_att_text(ncid, varid_date, "units", tunits);

    for (i = 0; i < obs->nobstypes; ++i)
        ncw_put_att_int(ncid, varid_type, obs->obstypes[i].name, 1, &i);

    for (i = 0; i < st_getsize(obs->products); ++i)
        ncw_put_att_int(ncid, varid_product, st_findstringbyindex(obs->products, i), 1, &i);

    for (i = 0; i < st_getsize(obs->instruments); ++i)
        ncw_put_att_int(ncid, varid_instrument, st_findstringbyindex(obs->instruments, i), 1, &i);

    for (i = 0; i < st_getsize(obs->datafiles); ++i) {
        char attname[NC_MAX_NAME];
        char* datafile = st_findstringbyindex(obs->datafiles, i);

        snprintf(attname, NC_MAX_NAME, "%d", i);
        ncw_put_att_text(ncid, varid_fid, attname, datafile);
    }

    if (obs->nccompression > 0)
        ncw_def_deflate(ncid, 0, 1, obs->nccompression);
    ncw_enddef(ncid);

    if (nobs == 0) {
        ncw_close(ncid);
        return;
    }

    id = malloc(nobs * sizeof(int));
    id_orig = malloc(nobs * sizeof(int));
    type = malloc(nobs * sizeof(short int));
    product = malloc(nobs * sizeof(short int));
    instrument = malloc(nobs * sizeof(short int));
    fid = malloc(nobs * sizeof(short int));
    batch = malloc(nobs * sizeof(int));
    value = malloc(nobs * sizeof(double));
    std = malloc(nobs * sizeof(double));
    lon = malloc(nobs * sizeof(double));
    lat = malloc(nobs * sizeof(double));
    depth = malloc(nobs * sizeof(double));
    model_depth = malloc(nobs * sizeof(double));
    fi = malloc(nobs * sizeof(double));
    fj = malloc(nobs * sizeof(double));
    fk = malloc(nobs * sizeof(double));
    date = malloc(nobs * sizeof(double));
    status = malloc(nobs * sizeof(int));
    aux = malloc(nobs * sizeof(int));

    for (i = 0, ii = 0; i < obs->nobs; ++i) {
        observation* m = &obs->data[i];

        if (!isfinite(m->value) || fabs(m->value) > FLT_MAX || !isfinite((float) m->value))
            enkf_quit("bad value");

        type[ii] = m->type;
        product[ii] = m->product;
        instrument[ii] = m->instrument;
        id[ii] = m->id;
        fid[ii] = m->fid;
        batch[ii] = m->batch;
        /*
         * id of the first ob contributed to this sob 
         */
        id_orig[ii] = m->id_orig;
        value[ii] = m->value;
        std[ii] = m->std;
        lon[ii] = m->lon;
        lat[ii] = m->lat;
        depth[ii] = m->depth;
        model_depth[ii] = m->model_depth;
        fi[ii] = m->fi;
        fj[ii] = m->fj;
        fk[ii] = m->fk;
        date[ii] = m->date;
        status[ii] = m->status;
        aux[ii] = m->aux;
        ii++;
    }
    assert(ii == nobs);

    ncw_put_var_int(ncid, varid_id, id);
    ncw_put_var_int(ncid, varid_idorig, id_orig);
    ncw_put_var_short(ncid, varid_type, type);
    ncw_put_var_short(ncid, varid_product, product);
    ncw_put_var_short(ncid, varid_instrument, instrument);
    ncw_put_var_short(ncid, varid_fid, fid);
    ncw_put_var_int(ncid, varid_batch, batch);
    ncw_put_var_double(ncid, varid_value, value);
    ncw_put_var_double(ncid, varid_std, std);
    ncw_put_var_double(ncid, varid_lon, lon);
    ncw_put_var_double(ncid, varid_lat, lat);
    ncw_put_var_double(ncid, varid_depth, depth);
    ncw_put_var_double(ncid, varid_mdepth, model_depth);
    ncw_put_var_double(ncid, varid_fi, fi);
    ncw_put_var_double(ncid, varid_fj, fj);
    ncw_put_var_double(ncid, varid_fk, fk);
    ncw_put_var_double(ncid, varid_date, date);
    ncw_put_var_int(ncid, varid_status, status);
    ncw_put_var_int(ncid, varid_aux, aux);

    ncw_close(ncid);
    free(type);
    free(product);
    free(instrument);
    free(id);
    free(id_orig);
    free(fid);
    free(batch);
    free(value);
    free(std);
    free(lon);
    free(lat);
    free(depth);
    free(model_depth);
    free(fi);
    free(fj);
    free(fk);
    free(date);
    free(status);
    free(aux);
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/**
 */
void obs_writeaux(observations* obs, char fname[])
{
    int ncid;
    int dimid_nobs;
    int varid_aux;
    size_t n;
    int i, ii;
    int* aux;

    if (rank != 0)
        return;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_varid(ncid, "aux", &varid_aux);
    ncw_inq_dimid(ncid, "nobs", &dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs, &n);
    assert((int) n == obs->nobs);

    aux = malloc(n * sizeof(int));
    for (i = 0, ii = 0; i < obs->nobs; ++i) {
        aux[ii] = obs->data[i].aux;
        ii++;
    }

    ncw_put_var_int(ncid, varid_aux, aux);
    ncw_close(ncid);
    free(aux);
}

/**
 */
int obs_modifiederrors_alreadywritten(observations* obs, char fname[])
{
    int ncid;
    int iswritten;

    ncw_open(fname, NC_NOWRITE, &ncid);
    iswritten = ncw_var_exists(ncid, "std_orig");
    ncw_close(ncid);

    return iswritten;
}

/**
 */
static int cmp_xyz(const void* p1, const void* p2)
{
    observation* o1 = (observation*) p1;
    observation* o2 = (observation*) p2;

    if (o1->lon > o2->lon)
        return 1;
    else if (o1->lon < o2->lon)
        return -1;
    else if (o1->lat > o2->lat)
        return 1;
    else if (o1->lat < o2->lat)
        return -1;
    else if (o1->depth > o2->depth)
        return 1;
    else if (o1->depth < o2->depth)
        return -1;
    else if (o1->instrument > o2->instrument)
        return 1;
    else if (o1->instrument < o2->instrument)
        return -1;
    return 0;
}

/**
 */
void obs_superob(observations* obs, __compar_d_fn_t cmp_obs, observations** sobs, int sobid, int do_thin)
{
    int i1 = 0, i2 = 0;
    int nsobs = 0;
    int nthinned = 0;
    observation* data = obs->data;
    observation* sdata = NULL;

    qsort_r(obs->data, obs->ngood, sizeof(observation), cmp_obs, obs);

    while (i2 < obs->ngood) {
        observation* so;
        observation* o;
        double lon_min, lon_max;
        int ii;
        double sum, sumsq, subvar;
        double var;

        /*
         * identify obs that will be combined into this superob 
         */
        while (i2 + 1 < obs->ngood && cmp_obs(&data[i1], &data[i2 + 1], obs) == 0)
            i2++;

        /*
         * thin observations with identical positions (these are supposedly
         * coming from high-frequency instruments)
         */
        if (do_thin) {
            int i11 = i1;
            int i22 = i1;

            qsort(&data[i1], i2 - i1 + 1, sizeof(observation), cmp_xyz);
            while (i22 <= i2) {
                while (i22 + 1 <= i2 && cmp_xyz(&data[i11], &data[i22 + 1]) == 0)
                    i22++;
                for (ii = i11 + 1; ii <= i22; ++ii) {
                    data[ii].status = STATUS_THINNED;
                    nthinned++;
                }
                i22++;
                i11 = i22;
            }
        }

        if (nsobs % NOBS_INC == 0)
            sdata = realloc(sdata, (nsobs + NOBS_INC) * sizeof(observation));

        if (sobid == nsobs) {
            int i;

            enkf_printf("    sob # %d is formed by the following obs:\n", sobid);
            for (i = i1; i <= i2; ++i) {
                enkf_printf("      ");
                obs_printob(obs, i);
            }
        }

        so = &sdata[nsobs];
        o = &data[i1];

        assert(o->status == STATUS_OK);

        if (!enkf_considersubgridvar)
            subvar = 0.0;
        else {
            /*
             * Calculate subgrid variance. For now assume equal weights.
             */
            sum = 0.0;
            sumsq = 0.0;
            for (ii = i1; ii <= i2; ++ii) {
                sum += data[ii].value;
                sumsq += data[ii].value * data[ii].value;
            }
            subvar = (i2 - i1 > 0) ? (sumsq - sum * sum / (double) (i2 - i1 + 1)) / (double) (i2 - i1) : 0.0;
        }

        /*
         * set the "aux" field in the original obs to the id of the superob
         * they are merged into
         */
        for (ii = i1; ii <= i2; ++ii)
            data[ii].aux = nsobs;

        so->type = o->type;
        so->product = o->product;
        so->instrument = o->instrument;
        so->id = nsobs;
        so->id_orig = o->id_orig;
        so->fid = o->fid;
        so->batch = o->batch;

        var = o->std * o->std;
        var = (subvar > var) ? subvar : var;
        so->std = 1.0 / var;
        so->value = o->value;
        so->lon = o->lon;
        so->lat = o->lat;
        so->depth = o->depth;
        so->model_depth = o->model_depth;
        so->fi = o->fi;
        so->fj = o->fj;
        so->fk = o->fk;
        so->date = o->date;
        so->status = o->status;
        so->aux = 1;

        lon_min = o->lon;
        lon_max = o->lon;
        for (ii = i1 + 1; ii <= i2; ++ii) {
            o = &data[ii];
            if (o->status == STATUS_THINNED)
                continue;
            if (so->product != o->product)
                so->product = -1;
            if (so->instrument != o->instrument)
                so->instrument = -1;
            if (so->fid != o->fid)
                so->fid = -1;
            if (so->batch != o->batch)
                so->batch = -1;

            var = o->std * o->std;
            if (subvar > var) {
                var = subvar;
                obs->obstypes[o->type].nsubgrid++;
            }
            so->value = so->value * so->std + o->value / var;
            so->lon = so->lon * so->std + o->lon / var;
            so->lat = so->lat * so->std + o->lat / var;
            so->depth = so->depth * so->std + o->depth / var;
            so->model_depth = so->model_depth * so->std + o->model_depth / var;
            so->fi = so->fi * so->std + o->fi / var;
            so->fj = so->fj * so->std + o->fj / var;
            so->fk = so->fk * so->std + o->fk / var;
            so->date = so->date * so->std + o->date / var;

            so->std += 1.0 / var;

            so->value /= so->std;
            so->lon /= so->std;
            so->lat /= so->std;
            so->depth /= so->std;
            so->model_depth /= so->std;
            so->fi /= so->std;
            so->fj /= so->std;
            so->fk /= so->std;
            so->date /= so->std;
            so->aux++;

            if (o->lon < lon_min)
                lon_min = o->lon;
            if (o->lon > lon_max)
                lon_max = o->lon;

            assert(o->status == STATUS_OK);
        }

        if (lon_max - lon_min > 180.0) {
            /*
             * (there is a possibility of merging observations separated by
             * more than 360 degrees in longitude) 
             */
            so->lon = 0.0;
            for (ii = i1; ii < i2; ++ii) {
                o = &data[ii];
                if (o->lon - lon_min > 180.0)
                    o->lon -= 360.0;
                var = o->std * o->std;
                var = (subvar > var) ? subvar : var;
                so->lon += o->lon / var;
            }
            so->lon *= so->std;
            if (so->lon < 0.0)
                so->lon += 360.0;
        }
        so->std = sqrt(1.0 / so->std);
        if (so->std < obs->obstypes[so->type].estdmin)
            so->std = obs->obstypes[so->type].estdmin;

        nsobs++;

        i1 = i2 + 1;
        i2 = i1;
    }
    if (nthinned > 0) {
        int i;

        for (i = 0; i < obs->nobs; ++i) {
            observation* m = &obs->data[i];
            obstype* ot = &obs->obstypes[m->type];

            if (m->status == STATUS_THINNED) {
                obs->nthinned++;
                ot->nthinned++;
            }
        }
        assert(nthinned == obs->nthinned);
        enkf_printf("    thinned %d observations\n", nthinned);
    }
    enkf_printf("    %d superobservations\n", nsobs);

    *sobs = obs_create_fromdata(obs, nsobs, sdata);
    if (sobid >= 0) {
        enkf_printf("    sob # %d info:\n", sobid);
        enkf_printf("      ");
        obs_printob(*sobs, sobid);
    }

    obs_calcstats(*sobs);
}

/**
 */
void obs_find_bytype(observations* obs, int type, int* nobs, int** obsids)
{
    int i;

    if (!enkf_fstatsonly)
        assert(obs->obstypes[type].nobs == obs->obstypes[type].ngood);

    *nobs = 0;
    if (obs->obstypes[type].nobs == 0) {
        *obsids = NULL;
        return;
    }
    *obsids = malloc(obs->obstypes[type].nobs * sizeof(int));
    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];

        if (o->type == type && o->status == STATUS_OK) {
            (*obsids)[*nobs] = i;
            (*nobs)++;
        }
    }
    if (*nobs == 0) {
        free(*obsids);
        *obsids = NULL;
    }
}

/**
 */
void obs_find_bytypeandtime(observations* obs, int type, int time, int* nobs, int** obsids)
{
    obstype* ot = &obs->obstypes[type];
    int i;

    if (!enkf_fstatsonly)
        assert(ot->nobs == ot->ngood);

    *nobs = 0;
    if (ot->ngood == 0) {
        *obsids = NULL;
        return;
    }
    *obsids = malloc(obs->obstypes[type].ngood * sizeof(int));
    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];

        if (o->type == type && o->status == STATUS_OK && get_tshift(o->date, ot->async_tstep, ot->async_centred) == time) {
            (*obsids)[*nobs] = i;
            (*nobs)++;
        }
    }
    if (*nobs == 0) {
        free(*obsids);
        *obsids = NULL;
    }
}

/**
 */
void obs_printob(observations* obs, int i)
{
    observation* o = &obs->data[i];

    enkf_printf("type = %s, product = %s, instrument = %s, datafile = %s, id = %d, original id = %d, batch = %d, value = %.3g, std = %.3g, ", obs->obstypes[o->type].name, st_findstringbyindex(obs->products, o->product), st_findstringbyindex(obs->instruments, o->instrument), st_findstringbyindex(obs->datafiles, o->fid), o->id, o->id_orig, (int) o->batch, o->value, o->std);
    enkf_printf("lon = %.3f, lat = %.3f, depth = %.1f, model_depth = %.1f, fi = %.3f, fj = %.3f, fk = %.3f, date = %.3g, status = %d\n", o->lon, o->lat, o->depth, o->model_depth, o->fi, o->fj, o->fk, o->date, o->status);
}

#if defined(ENKF_CALC)
/**
 */
static double distance(double xyz1[3], double xyz2[3])
{
    return sqrt((xyz1[0] - xyz2[0]) * (xyz1[0] - xyz2[0]) + (xyz1[1] - xyz2[1]) * (xyz1[1] - xyz2[1]) + (xyz1[2] - xyz2[2]) * (xyz1[2] - xyz2[2]));
}

/**
 */
void obs_createkdtrees(observations* obs, model* m)
{
    int otid;

    if (obs->loctrees == NULL)
        obs->loctrees = calloc(obs->nobstypes, sizeof(kdtree*));
    if (obs->obsids == NULL)
        obs->obsids = calloc(obs->nobstypes, sizeof(int*));

    for (otid = 0; otid < obs->nobstypes; ++otid) {
        obstype* ot = &obs->obstypes[otid];
        grid* g = model_getgridbyid(m, ot->gridid);
        kdtree** tree = &obs->loctrees[otid];
        int nobs = 0;
        int* obsids = NULL;

#if defined(OBS_SHUFFLE)
        int* ids = NULL;
#endif
        int i;

        obs_find_bytype(obs, otid, &nobs, &obsids);
        if (nobs == 0)
            continue;
        ot->nobs = nobs;
        if (obs->obsids[otid] != NULL)
            free(obs->obsids[otid]);
        obs->obsids[otid] = obsids;

        if (*tree != NULL)
            kd_destroy(*tree);

#if defined(OBS_SHUFFLE)
        ids = malloc(nobs * sizeof(int));
        for (i = 0; i < nobs; ++i)
            ids[i] = i;
        shuffle(nobs, ids);
#endif

        *tree = kd_create(3);
        for (i = 0; i < nobs; ++i) {
#if defined(OBS_SHUFFLE)
            int id = ids[i];
            observation* o = &obs->data[obsids[id]];
#else
            observation* o = &obs->data[obsids[i]];
#endif
            double ll[2] = { o->lon, o->lat };
            double xyz[3];

            grid_tocartesian(g, ll, xyz);
#if defined(OBS_SHUFFLE)
            kd_insertnode(*tree, xyz, id);
#else
            kd_insertnode(*tree, xyz, i);
#endif
        }

#if defined(OBS_SHUFFLE)
        free(ids);
#endif
    }
}

/**
 */
void obs_destroykdtrees(observations* obs)
{
    int otid;

    for (otid = 0; otid < obs->nobstypes; ++otid) {
        kdtree** tree = &obs->loctrees[otid];

        if (*tree != NULL)
            kd_destroy(*tree);
        *tree = NULL;
    }
}

/**
 */
#if defined(MINIMISE_ALLOC)
void obs_findlocal(observations* obs, model* m, grid* g, int icoord, int jcoord, int* n, int** ids, double** lcoeffs, int* ploc_allocated)
#else
void obs_findlocal(observations* obs, model* m, grid* g, int icoord, int jcoord, int* n, int** ids, double** lcoeffs)
#endif
{
    double ll[2];
    double xyz[3];
    int otid;
    int i, ntot, ngood;

    grid_ij2xy(g, icoord, jcoord, &ll[0], &ll[1]);
    grid_tocartesian(g, ll, xyz);

    if (obs->nobstypes == 0)
        return;
    if (obs->loctrees == NULL)
        obs_createkdtrees(obs, m);

    for (otid = 0, i = 0; otid < obs->nobstypes; ++otid) {
        obstype* ot = &obs->obstypes[otid];
        kdtree* tree = obs->loctrees[otid];
        int* obsids = obs->obsids[otid];
        kdset* set = NULL;
        double dist;
        size_t id;
        int iloc;

        if (ot->nobs == 0)
            continue;

        if (ot->ndomains > 0) {
            /*
             * (if ot->ndomains = 0 then observations of this type are visible
             * from all grids)
             */
            int d;

            for (d = 0; d < ot->ndomains; ++d)
                if (strcasecmp(grid_getdomainname(g), ot->domainnames[d]) == 0)
                    break;
            if (d == ot->ndomains)
                continue;
        }

        set = kd_findnodeswithinrange(tree, xyz, obstype_getmaxlocrad(ot), 1);
        for (iloc = 0; iloc < ot->nlobsmax && (id = kdset_readnext(set, &dist)) != SIZE_MAX; ++i, ++iloc) {
            int id_orig = kd_getnodeorigid(tree, id);

#if defined(MINIMISE_ALLOC)
            if (ploc_allocated != NULL) {
                if (i >= *ploc_allocated) {
                    *ploc_allocated += KD_INC;
                    *ids = realloc(*ids, *ploc_allocated * sizeof(int));
                    *lcoeffs = realloc(*lcoeffs, *ploc_allocated * sizeof(double));
                }
            } else {
                if (i % KD_INC == 0) {
                    *ids = realloc(*ids, (i + KD_INC) * sizeof(int));
                    *lcoeffs = realloc(*lcoeffs, (i + KD_INC) * sizeof(double));
                }
            }
#else
            if (i % KD_INC == 0) {
                *ids = realloc(*ids, (i + KD_INC) * sizeof(int));
                *lcoeffs = realloc(*lcoeffs, (i + KD_INC) * sizeof(double));
            }
#endif
            (*ids)[i] = obsids[id_orig];
        }
        kdset_free(set);
    }

    /*
     * compact the result, calculate taper coefficients
     */
    ntot = i;
    for (i = 0, ngood = 0; i < ntot; ++i) {
        int id = (*ids)[i];
        observation* o = &obs->data[id];
        double ll2[2] = { o->lon, o->lat };
        double xyz2[3];

        if (o->status != STATUS_OK)
            continue;

        grid_tocartesian(g, ll2, xyz2);
        (*ids)[ngood] = id;
        (*lcoeffs)[ngood] = obstype_calclcoeff(&obs->obstypes[o->type], distance(xyz, xyz2));
        ngood++;
    }
    *n = ngood;
#if defined (MINIMISE_ALLOC)
    if (ploc_allocated == NULL && ngood == 0 && *ids != NULL) {
#else
    if (ngood == 0 && *ids != NULL) {
#endif
        free(*ids);
        free(*lcoeffs);
    }
}
#endif
