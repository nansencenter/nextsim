/******************************************************************************
 *
 * File:        reader_cs2smos.c        
 *
 * Created:     6/2020
 *
 * Author:      Sukun Cheng
 *              NERSC
 *
 * Description: Contains reader_cs2smos_standard() for preprocessed
 *              data from CS2-SMOS.
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


/**
 * read CS2_SMOS_v2.2 merged sea ice thickness, two dimensions, 
   W_XX-ESA,SMOS_CS2,NH_25KM_EASE2_20181111_20181117_r_v202_01_l4sit.nc << v202
*/
void reader_cs2smos_standard(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    int ncid;
    int x_dimid, y_dimid, t_dimid, lon_dimid, lat_dimid;
    size_t nobslen, xlen, ylen, tlen, lonlen, latlen;
    int varid_lon, varid_lat, varid_sit, varid_error, varid_sic, varid_time;
    int   sit_fill_value,  estd_fill_value,  sic_fill_value;
    float sit_scale_factor,estd_scale_factor,sic_scale_factor;
    float** lon = NULL;
    float** lat = NULL;
    float* coast_lon = NULL;
    float* coast_lat = NULL;
    int*** sit = NULL;
    int*** error_std = NULL;
    int*** sic = NULL;
    double* time = NULL;
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
    lat = alloc2d(ylen, xlen, sizeof(float));  // use 2d array. sit[i*j+j]] in standard2 is a bug to be fixed 
    ncw_get_var_float(ncid, varid_lat, lat[0]);
    ncw_inq_varid(ncid, "lon", &varid_lon);
    lon = alloc2d(ylen, xlen, sizeof(float));
    ncw_get_var_float(ncid, varid_lon, lon[0]);
    error_std = alloc3d(tlen, ylen, xlen, sizeof(int));
    sit = alloc3d(tlen, ylen, xlen, sizeof(int));
    sic = alloc3d(tlen, ylen, xlen, sizeof(int));
    
    if (strcmp(meta->type,"sea_ice_thickness") == 0) {  // type in obstypes.prm
        ncw_inq_varid(ncid, "analysis_sea_ice_thickness", &varid_sit);     
        ncw_get_var_int(ncid, varid_sit, sit[0][0]);
        if (ncw_att_exists(ncid, varid_sit, "_FillValue"))
            ncw_get_att_int(ncid, varid_sit, "_FillValue", &sit_fill_value);
        if (ncw_att_exists(ncid, varid_sit, "scale_factor"))
            ncw_get_att_float(ncid, varid_sit, "scale_factor", &sit_scale_factor);

        ncw_inq_varid(ncid, "analysis_sea_ice_thickness_unc", &varid_error);    
        ncw_get_var_int(ncid, varid_error, error_std[0][0]);
        if (ncw_att_exists(ncid, varid_error, "_FillValue"))
            ncw_get_att_int(ncid, varid_error, "_FillValue", &estd_fill_value);
        if (ncw_att_exists(ncid, varid_error, "scale_factor"))
            ncw_get_att_float(ncid, varid_error, "scale_factor", &estd_scale_factor);
    }
    else if (strcmp(meta->type,"sea_ice_concentration") == 0) {
        ncw_inq_varid(ncid, "sea_ice_concentration", &varid_sic);
        ncw_get_var_int(ncid, varid_sic, sic[0][0]);
        if (ncw_att_exists(ncid, varid_sic, "_FillValue"))
            ncw_get_att_int(ncid, varid_sic, "_FillValue", &sic_fill_value);
        if (ncw_att_exists(ncid, varid_sic, "scale_factor"))
            ncw_get_att_float(ncid, varid_sic, "scale_factor", &sic_scale_factor);        
    }
    ncw_inq_varid(ncid, "time", &varid_time);
    time = malloc(tlen * sizeof(double));
    ncw_get_var_double(ncid, varid_time, time);
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

    // load coast data from 'reference_grid_coast.nc
    fname='reference_grid_coast.nc';
    ncw_open(fname, NC_NOWRITE, &ncid);
    // longitude
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "lon")) ? "lon" : "length", &lon_dimid);
    ncw_inq_dimlen(ncid, lon_dimid, &lonlen);
    ncw_inq_varid(ncid, "lon", &varid_lon);
    coast_lon = malloc(lonlen * sizeof(float));
    // latitude
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "lat")) ? "lat" : "length", &lat_dimid);
    ncw_inq_dimlen(ncid, lat_dimid, &latlen);
    ncw_inq_varid(ncid, "lat", &varid_lat);
    coast_lat = malloc(latlen * sizeof(float));
    ncw_close(ncid);
    
    //
    for (int it= 0; it< (int) tlen; ++it) {
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
                if (strcmp(meta->type,"sea_ice_thickness") == 0) {
                    o->value = (double) (sit[it][i][j]*sit_scale_factor);
                    // o->std   = (double) (error_std[it][i][j]*estd_scale_factor);
                    o->std   = (double) (fmin(0.5, 0.1 + 0.15*sit[it][i][j])*estd_scale_factor);
                }
                else if (strcmp(meta->type,"sea_ice_concentration") == 0) {
                    o->value = (double) (sic[it][i][j]*sic_scale_factor*0.01);
                    o->std   = 0.01+pow(0.5-fabs(0.5-o->value),2.0);
                }
                o->lon = lon[i][j];
                o->lat = lat[i][j];
                o->depth = 0.0;
                o->fk = (double) ksurf;
                o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
                if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                    continue;
                // check distance between the observation the closest coast
                o->status = obs_distance2coast(g, coast_lon, coast_lat, o->lon, o->lat); 
                if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                    continue;
                o->model_depth = NAN;   // set in obs_add()
                o->date = time[tlen] * tunits_multiple + tunits_offset;
                o->aux = -1;

                obs->nobs++;
            }
        }
    }

    free(lon);
    free(lat);
    if (strcmp(meta->type,"sea_ice_thickness") == 0) {
        free(sit);
        free(error_std);
    }
    else if (strcmp(meta->type,"sea_ice_concentration") == 0) 
        free(sic);
    
    free(time);
}

/**
 * @brief      filter out observations near coast (min_distance2coast_km)
 * 
 * @param g    grid object, containing coast coordinates.
 * @param lon  longtitude of observation
 * @param lat  latitude of observation
 * @return int return status
 */
int obs_distance2coast(grid* g, float* coast_lon, float* coast_lat, float lon, float lat){ 
    int min_distance2coast_km = 30;
    double ll[2] = { lon, lat }; //observation coordinates
    double xyz1[3];
    grid_tocartesian(g, ll, xyz1);
    // distance filter
    for (int j = 0; j < sizeof(coast_lon); ++j) {
        double ll2[2] = { coast_lon[j], coast_lat[j] };
        double xyz2[3];    
        grid_tocartesian(g, ll2, xyz2);
        // reject observation, borrow STATUS_OUTSIDEGRID.
        double distance = sqrt((xyz1[0] - xyz2[0]) * (xyz1[0] - xyz2[0]) + (xyz1[1] -  xyz2[1]) * (xyz1[1] - xyz2[1]) + (xyz1[2] - xyz2[2]) * (xyz1[2] - xyz2[2]));
        if (distance<=min_distance2coast_km) 
            return STATUS_OUTSIDEGRID; 
    }
    // accept the observation.
    return 0;  
}
