/******************************************************************************
 *
 * File:        reader_smos.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Contains 2 readers reader_smos_standard() for preprocessed
 *              data from SMOS.
 *              Parameters:
 *               
 *
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
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"


/** read one dimensional data
 */
void reader_smos_standard(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    int ncid;
    int dimid_nobs;
    size_t nobs_local;
    int varid_lon, varid_lat, varid_sit, varid_error, varid_time;
    int sit_fill_value;
    int estd_fill_value;
    double* lon = NULL;
    double* lat = NULL;
    double* sit = NULL;
    double* error_std = NULL;
    double* time = NULL;
    int year, month, day;
    char tunits[MAXSTRLEN];
    size_t tunits_len;
    double tunits_multiple, tunits_offset;
    char* basename;
    int i;

    basename = strrchr(fname, '/');
    if (basename == NULL)
        basename = fname;
    else
        basename += 1;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "nobs")) ? "nobs" : "length", &dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs, &nobs_local);
    enkf_printf("        nobs = %u\n", (unsigned int) nobs_local);

    if (nobs_local == 0) {
        ncw_close(ncid);
        return;
    }

    ncw_inq_varid(ncid, "latitude", &varid_lat);
    lat = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "longitude", &varid_lon);
    lon = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_lon, lon);

    ncw_inq_varid(ncid, "sea_ice_thickness", &varid_sit);
    sit = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_sit, sit);
    ncw_get_att_int(ncid, varid_sit, "missing_value", &sit_fill_value);

    ncw_inq_varid(ncid, "ice_thickness_uncertainty", &varid_error);
    error_std = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_error, error_std);
    ncw_get_att_int(ncid, varid_sit, "missing_value", &estd_fill_value);

    ncw_inq_varid(ncid, "time", &varid_time);
    time = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_time, time);
    ncw_inq_attlen(ncid, varid_time, "units", &tunits_len);
    ncw_get_att_text(ncid, varid_time, "units", tunits);
    basename[37] = 0;
    if (!str2int(&basename[35], &day))
        enkf_quit("SMOS reader: could not convert file name \"%s\" to date", fname);
    basename[11] = 0;
    if (!str2int(&basename[35], &month))
        enkf_quit("SMOS reader: could not convert file name \"%s\" to date", fname);
    basename[9] = 0;
    if (!str2int(&basename[35], &year))
        enkf_quit("SMOS reader: could not convert file name \"%s\" to date", fname);
    snprintf(&tunits[tunits_len], MAXSTRLEN - tunits_len, " since %4d-%02d-%02d", year, month, day);

    ncw_close(ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    for (i = 0; i < (int) nobs_local; ++i) {
        observation* o;

        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        o->instrument = st_add_ifabsent(obs->instruments, "SMOS", -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = 0;
        o->value = sit[i];
        o->std = error_std[i];
        o->lon = lon[i];
        o->lat = lat[i];
        o->depth = 0.0;
        o->fk = (double) ksurf;
        o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
//        enkf_printf("   lon = %f lat = %f\n", o->lon, o->lat);
//        enkf_printf("   sit = %f\n", o->value);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
           continue;
        o->model_depth = NAN;   /* set in obs_add() */
        o->date = time[i] * tunits_multiple + tunits_offset;
        o->aux = -1;

        obs->nobs++;
    }

    free(lon);
    free(lat);
    free(sit);
    free(error_std);
    free(time);
}
/** read two dimensional data
 */
void reader_smos_standard2(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    int ncid;
    int x_dimid, y_dimid, t_dimid;
    size_t nobslen, xlen, ylen, tlen;
    int varid_lon, varid_lat, varid_sit, varid_error, varid_time;
    int sit_fill_value;
    int estd_fill_value;
    float* lon = NULL;
    float* lat = NULL;
    float* sit = NULL;
    float* error_std = NULL;
    float* time = NULL;
    int year, month, day;
    char tunits[MAXSTRLEN];
    size_t tunits_len;
    double tunits_multiple, tunits_offset;
    char* basename;

    basename = strrchr(fname, '/');
    if (basename == NULL)
        basename = fname;
    else
        basename += 1;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "x")) ? "x" : "length", &x_dimid);
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "y")) ? "y" : "length", &y_dimid);
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "time")) ? "time" : "length", &t_dimid);
    ncw_inq_dimlen(ncid, x_dimid, &xlen);
    ncw_inq_dimlen(ncid, y_dimid, &ylen);
    ncw_inq_dimlen(ncid, t_dimid, &tlen);
    nobslen = xlen * ylen;
    enkf_printf("        x = %u\n y = %u\n t = %u\n nobs = %u\n", xlen, ylen, tlen, nobslen);

    ncw_inq_varid(ncid, "latitude", &varid_lat);
    lat = malloc(ylen*xlen*sizeof(float));
    ncw_get_var_float(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "longitude", &varid_lon);
    lon = malloc(ylen*xlen*sizeof(float));
    ncw_get_var_float(ncid, varid_lon, lon);

    ncw_inq_varid(ncid, "sea_ice_thickness", &varid_sit);
    sit = malloc(tlen*ylen*xlen*sizeof(float));
    ncw_get_var_float(ncid, varid_sit, sit);
    ncw_get_att_int(ncid, varid_sit, "missing_value", &sit_fill_value);

    ncw_inq_varid(ncid, "ice_thickness_uncertainty", &varid_error);
    error_std = malloc(tlen*ylen*xlen*sizeof(float));
    ncw_get_var_float(ncid, varid_error, error_std);
    ncw_get_att_int(ncid, varid_sit, "missing_value", &estd_fill_value);

    ncw_inq_varid(ncid, "time", &varid_time);
    time = malloc(tlen * sizeof(float));
    ncw_get_var_float(ncid, varid_time, time);
    ncw_inq_attlen(ncid, varid_time, "units", &tunits_len);
    ncw_get_att_text(ncid, varid_time, "units", tunits);
    basename[37] = 0;
    if (!str2int(&basename[35], &day))
        enkf_quit("SMOS reader: could not convert file name \"%s\" to date", fname);
    basename[11] = 0;
    if (!str2int(&basename[35], &month))
        enkf_quit("SMOS reader: could not convert file name \"%s\" to date", fname);
    basename[9] = 0;
    if (!str2int(&basename[35], &year))
        enkf_quit("SMOS reader: could not convert file name \"%s\" to date", fname);
    snprintf(&tunits[tunits_len], MAXSTRLEN - tunits_len, " since %4d-%02d-%02d", year, month, day);

    ncw_close(ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    for (int i = 0; i < (int) ylen; ++i) {
        for (int j = 0; j < (int) xlen; ++j) {

            observation* o;

            obs_checkalloc(obs);
            o = &obs->data[obs->nobs];

            o->product = st_findindexbystring(obs->products, meta->product);
            assert(o->product >= 0);
            o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
            o->instrument = st_add_ifabsent(obs->instruments, "SMOS", -1);
            o->id = obs->nobs;
            o->fid = fid;
            o->batch = 0;
            o->value = sit[i*j+j];
//            o->std = error_std[i*j+j]; # there are zeros in the smos data for std. Not working for EnKF
            o->std = .02;
            o->lon = lon[i*j+j];
            o->lat = lat[i*j+j];
            o->depth = 0.0;
            o->fk = (double) ksurf;
            o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
//            enkf_printf("   lon = %f lat = %f\n val= %f\n std = %f\n", o->lon, o->lat, o->value, o->std);
//            enkf_printf("   sit = %f\n", o->value);
            if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
               continue;
            o->model_depth = NAN;   // set in obs_add()
            o->date = time[tlen] * tunits_multiple + tunits_offset;
            o->aux = -1;

            obs->nobs++;
        }
    }

    free(lon);
    free(lat);
    free(sit);
    free(error_std);
    free(time);
}

/**
 * read cs2-smos merged sea ice thickness, two dimensions
*/
void reader_smos_standard3(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    int ncid;
    int x_dimid, y_dimid, t_dimid;
    size_t nobslen, xlen, ylen, tlen;
    int varid_lon, varid_lat, varid_sit, varid_error, varid_time;
    int sit_fill_value;
    int estd_fill_value;
    float* lon = NULL;
    float* lat = NULL;
    float* sit = NULL;
    float* error_std = NULL;
    float* time = NULL;
    int year, month, day;
    char tunits[MAXSTRLEN];
    size_t tunits_len;
    double tunits_multiple, tunits_offset;
    char* basename;

    basename = strrchr(fname, '/');
    if (basename == NULL)
        basename = fname;
    else
        basename += 1;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "xc")) ? "xc" : "length", &x_dimid);
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "yc")) ? "yc" : "length", &y_dimid);
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "time")) ? "time" : "length", &t_dimid);
    ncw_inq_dimlen(ncid, x_dimid, &xlen);
    ncw_inq_dimlen(ncid, y_dimid, &ylen);
    ncw_inq_dimlen(ncid, t_dimid, &tlen);
    nobslen = xlen * ylen;
    enkf_printf("        x = %u\n y = %u\n t = %u\n nobs = %u\n", xlen, ylen, tlen, nobslen);

    ncw_inq_varid(ncid, "lat", &varid_lat);
    lat = malloc(ylen*xlen*sizeof(float));
    ncw_get_var_float(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "lon", &varid_lon);
    lon = malloc(ylen*xlen*sizeof(float));
    ncw_get_var_float(ncid, varid_lon, lon);

    ncw_inq_varid(ncid, "analysis_sea_ice_thickness", &varid_sit);
    sit = malloc(tlen*ylen*xlen*sizeof(float));
    ncw_get_var_float(ncid, varid_sit, sit);
    if (ncw_att_exists(ncid, varid_sit, "_FillValue"))
        ncw_get_att_int(ncid, varid_sit, "_FillValue", &sit_fill_value);
    else if (ncw_att_exists(ncid, varid_sit, "missing_value"))
        ncw_get_att_int(ncid, varid_sit, "missing_value", &sit_fill_value);

    ncw_inq_varid(ncid, "analysis_sea_ice_thickness_unc", &varid_error);
    error_std = malloc(tlen*ylen*xlen*sizeof(float));
    ncw_get_var_float(ncid, varid_error, error_std);
    if (ncw_att_exists(ncid, varid_error, "_FillValue"))
        ncw_get_att_int(ncid, varid_error, "_FillValue", &estd_fill_value);
    else if (ncw_att_exists(ncid, varid_error, "missing_value"))
        ncw_get_att_int(ncid, varid_error, "missing_value", &estd_fill_value);

    ncw_inq_varid(ncid, "time", &varid_time);
    time = malloc(tlen * sizeof(float));
    ncw_get_var_float(ncid, varid_time, time);
    ncw_inq_attlen(ncid, varid_time, "units", &tunits_len);
    ncw_get_att_text(ncid, varid_time, "units", tunits);
    basename[37] = 0;
    if (!str2int(&basename[35], &day))
        enkf_quit("CS2SMOS reader: could not convert file name \"%s\" to date", fname);
    basename[11] = 0;
    if (!str2int(&basename[35], &month))
        enkf_quit("CS2SMOS reader: could not convert file name \"%s\" to date", fname);
    basename[9] = 0;
    if (!str2int(&basename[35], &year))
        enkf_quit("SMOS reader: could not convert file name \"%s\" to date", fname);
    snprintf(&tunits[tunits_len], MAXSTRLEN - tunits_len, " since %4d-%02d-%02d", year, month, day);

    ncw_close(ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    for (int i = 0; i < (int) ylen; ++i) {
        for (int j = 0; j < (int) xlen; ++j) {

            observation* o;

            obs_checkalloc(obs);
            o = &obs->data[obs->nobs];

            o->product = st_findindexbystring(obs->products, meta->product);
            assert(o->product >= 0);
            o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
            o->instrument = st_add_ifabsent(obs->instruments, "SMOS", -1);
            o->id = obs->nobs;
            o->fid = fid;
            o->batch = 0;
            o->value = sit[i*j+j];
//            o->std = error_std[i*j+j]; # there are zeros in the smos data for std. Not working for EnKF
            o->std = .02;
            o->lon = lon[i*j+j];
            o->lat = lat[i*j+j];
            o->depth = 0.0;
            o->fk = (double) ksurf;
            o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
//            enkf_printf("   lon = %f lat = %f\n val= %f\n std = %f\n", o->lon, o->lat, o->value, o->std);
//            enkf_printf("   sit = %f\n", o->value);
            if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
               continue;
            o->model_depth = NAN;   // set in obs_add()
            o->date = time[tlen] * tunits_multiple + tunits_offset;
            o->aux = -1;

            obs->nobs++;
        }
    }

    free(lon);
    free(lat);
    free(sit);
    free(error_std);
    free(time);
}
