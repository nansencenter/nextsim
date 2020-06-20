/******************************************************************************
 *
 * File:        reader_rads.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Contains 2 readers.
 *              - reader_rads_standard() -- for files of the form
 *                ??_yyyymmdd.nc. These files are assumed to have variable
 *                "time".
 *                Parameters:
 *                  MINDEPTH -- minimal allowed depth in metres.
 *              - reader_rads_standard2() -- for files of the form
 *                y<yyyy>/m<mm>/??_d<dd>.nc with no "time" variable.
 *                Parameters:
 *                  MINDEPTH -- minimal allowed depth in metres.
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "grid.h"
#include "model.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

/** For files of the form ??_yyyymmdd.nc. They are assumed to have "time" 
 * variable.
 */
void reader_rads_standard(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    char* addname = NULL;
    int ncid;
    int dimid_nobs;
    size_t nobs_local;
    int varid_lon, varid_lat, varid_pass, varid_sla, varid_time, varid_add;
    double* lon;
    double* lat;
    int* pass;
    double* sla;
    double* time;
    double error_std;
    double* add = NULL;
    size_t tunits_len;
    char* tunits;
    double tunits_multiple, tunits_offset;
    char* basename;
    char instname[3];
    int i;

    for (i = 0; i < meta->npars; ++i) {
        /*
         * (MINDEPTH is handled in obs_add() )
         */
        if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
            double mindepth;

            if (!str2double(meta->pars[i].value, &mindepth))
                enkf_quit("%s: can not convert MINDEPTH = \"%s\" to double\n", meta->prmfname, meta->pars[i].value);
            enkf_printf("        MINDEPTH = %.0f\n", mindepth);
            continue;
        } else if (strcasecmp(meta->pars[i].name, "ADD") == 0) {
            addname = meta->pars[i].value;
            enkf_printf("        ADDING \"%s\"\n", addname);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_inq_dimid(ncid, "nobs", &dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs, &nobs_local);
    enkf_printf("        nobs = %u\n", (unsigned int) nobs_local);

    if (nobs_local == 0) {
        ncw_close(ncid);
        return;
    }

    ncw_inq_varid(ncid, "lon", &varid_lon);
    lon = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_lon, lon);

    ncw_inq_varid(ncid, "lat", &varid_lat);
    lat = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "pass", &varid_pass);
    pass = malloc(nobs_local * sizeof(int));
    ncw_get_var_int(ncid, varid_pass, pass);

    ncw_inq_varid(ncid, "sla", &varid_sla);
    sla = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_sla, sla);
    ncw_get_att_double(ncid, varid_sla, "error_std", &error_std);
    enkf_printf("        error_std = %3g\n", error_std);

    if (addname != NULL) {
        ncw_inq_varid(ncid, addname, &varid_add);
        add = malloc(nobs_local * sizeof(double));
        ncw_get_var_double(ncid, varid_add, add);
    }

    ncw_inq_varid(ncid, "time", &varid_time);
    time = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_time, time);
    ncw_inq_attlen(ncid, varid_time, "units", &tunits_len);
    tunits = calloc(tunits_len + 1, 1);
    ncw_get_att_text(ncid, varid_time, "units", tunits);

    ncw_close(ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    basename = strrchr(fname, '/');
    if (basename == NULL)
        basename = fname;
    else
        basename += 1;
    strncpy(instname, basename, 2);
    instname[2] = 0;

    for (i = 0; i < (int) nobs_local; ++i) {
        observation* o;

        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        o->instrument = st_add_ifabsent(obs->instruments, instname, -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = pass[i];
        o->value = sla[i];
        if (add != NULL)
            o->value += add[i];
        o->std = error_std;
        o->lon = lon[i];
        o->lat = lat[i];
        o->depth = 0.0;
        o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        o->fk = (double) ksurf;
        o->model_depth = NAN;   /* set in obs_add() */
        o->date = time[i] * tunits_multiple + tunits_offset;

        o->aux = -1;

        obs->nobs++;
    }

    if (add != NULL)
        free(add);
    free(lon);
    free(lat);
    free(pass);
    free(sla);
    free(tunits);
    free(time);
}

/** For files of the form y<yyyy>/m<mm>/??_d<dd>.nc with no "time" variable.
 */
void reader_rads_standard2(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    char* addname = NULL;
    int ncid;
    int dimid_nobs;
    size_t nobs_local;
    int varid_lon, varid_lat, varid_pass, varid_sla, varid_flag, varid_add;
    double* lon;
    double* lat;
    int* pass;
    double* sla;
    int* flag;
    double error_std;
    double* add = NULL;
    char buf[MAXSTRLEN];
    int len;
    int year, month, day;
    double tunits_multiple, tunits_offset;
    char* basename;
    char instname[3];
    int i;

    for (i = 0; i < meta->npars; ++i) {
        /*
         * (MINDEPTH is handled in obs_add() )
         */
        if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
            double mindepth;

            if (!str2double(meta->pars[i].value, &mindepth))
                enkf_quit("observation prm file: can not convert MINDEPTH = \"%s\" to double\n", meta->pars[i].value);
            enkf_printf("        MINDEPTH = %.0f\n", mindepth);
            continue;
        } else if (strcasecmp(meta->pars[i].name, "ADD") == 0) {
            addname = meta->pars[i].value;
            enkf_printf("        ADDING \"%s\"\n", addname);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(ncid, "nobs", &dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs, &nobs_local);
    enkf_printf("        nobs = %u\n", (unsigned int) nobs_local);

    if (nobs_local == 0) {
        ncw_close(ncid);
        return;
    }

    ncw_inq_varid(ncid, "lon", &varid_lon);
    lon = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_lon, lon);

    ncw_inq_varid(ncid, "lat", &varid_lat);
    lat = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "pass", &varid_pass);
    pass = malloc(nobs_local * sizeof(int));
    ncw_get_var_int(ncid, varid_pass, pass);

    ncw_inq_varid(ncid, "sla", &varid_sla);
    sla = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_sla, sla);
    ncw_get_att_double(ncid, varid_sla, "error_std", &error_std);
    enkf_printf("        error_std = %3g\n", error_std);

    if (addname != NULL) {
        ncw_inq_varid(ncid, addname, &varid_add);
        add = malloc(nobs_local * sizeof(double));
        ncw_get_var_double(ncid, varid_add, add);
    }

    ncw_inq_varid(ncid, "local_flag", &varid_flag);
    flag = malloc(nobs_local * sizeof(int));
    ncw_get_var_int(ncid, varid_flag, flag);

    ncw_close(ncid);

    strcpy(buf, fname);
    len = strlen(buf);
    buf[len - 3] = 0;           /* .nc */
    if (!str2int(&buf[len - 5], &day))
        enkf_quit("RADS reader: could not convert file name \"%s\" to date", fname);
    buf[len - 10] = 0;
    if (!str2int(&buf[len - 12], &month))
        enkf_quit("RADS reader: could not convert file name \"%s\" to date", fname);
    buf[len - 14] = 0;
    if (!str2int(&buf[len - 18], &year))
        enkf_quit("RADS reader: could not convert file name \"%s\" to date", fname);
    snprintf(buf, MAXSTRLEN, "days since %4d-%02d-%02d", year, month, day);

    tunits_convert(buf, &tunits_multiple, &tunits_offset);

    basename = strrchr(fname, '/');
    if (basename == NULL)
        basename = fname;
    else
        basename += 1;
    strncpy(instname, basename, 2);
    instname[2] = 0;

    for (i = 0; i < (int) nobs_local; ++i) {
        observation* o;
        obstype* ot;

        if (flag[i] != 0)
            continue;

        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        ot = &obs->obstypes[o->type];
        o->instrument = st_add_ifabsent(obs->instruments, instname, -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = pass[i];
        o->value = sla[i];
        if (add != NULL)
            o->value += add[i];
        o->std = error_std;
        o->lon = lon[i];
        o->lat = lat[i];
        o->depth = 0.0;
        o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        o->model_depth = NAN;   /* set in obs_add() */
        o->fk = (double) ksurf;
        o->date = tunits_offset + 0.5;
        if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax))
            o->status = STATUS_OUTSIDEOBSDOMAIN;

        o->aux = -1;

        obs->nobs++;
    }

    if (add != NULL)
        free(add);
    free(lon);
    free(lat);
    free(pass);
    free(sla);
    free(flag);
}
