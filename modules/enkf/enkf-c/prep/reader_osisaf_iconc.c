/******************************************************************************
 *
 * File:        reader_osisaf_iconc.c        
 *
 * Created:     6/2020
 *
 * Author:      Sukun Cheng
 *              NERSC
 *
 * Description: read OSISAF ice concentration, product 401 (Global Sea Ice Concentration (SSMIS)), daily, 10km
 * http://www.osi-saf.org/?q=content/global-sea-ice-concentration-ssmis
 *               
 *
 *   Observations are loaded and selected within a domain defined by the reference grid: reference_grid.nc created from OCRA grid, used in NEMO. 
 *   The grid use variable mask to indicate land and ocean.
 *   mask is modified that the ocean area is limited with the nextsim domain by its output prior.nc sit, that mask(isnan(sit))=0
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


void reader_osisaf_iconc_standard(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    int ncid;
    int x_dimid, y_dimid, t_dimid, lon_dimid, lat_dimid;
    size_t nobslen, xlen, ylen, tlen, coast_len;
    int varid_lon, varid_lat, varid_sic, varid_time;
    int sic_fill_value;
    float  sic_scale_factor;
    float** lon = NULL;
    float** lat = NULL;
    float* coast_lon = NULL;
    float* coast_lat = NULL;
    int*** sic = NULL;
    double* time = NULL;
    int year, month, day;
    char* filename;
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
    enkf_printf(" reader_osisaf_iconc.c:\n");
    printf("         x = %lu\n         y = %lu\n         t = %lu\n         nobs =x*y= %lu\n", xlen, ylen, tlen, nobslen);
    ncw_inq_varid(ncid, "lat", &varid_lat);
    lat = alloc2d(ylen, xlen, sizeof(float));  
    ncw_get_var_float(ncid, varid_lat, lat[0]);
    ncw_inq_varid(ncid, "lon", &varid_lon);
    lon = alloc2d(ylen, xlen, sizeof(float));
    ncw_get_var_float(ncid, varid_lon, lon[0]);
    sic = alloc3d(tlen, ylen, xlen, sizeof(int));
    
    if (strcmp(meta->type,"SIC") == 0) {
        ncw_inq_varid(ncid, "ice_conc", &varid_sic);
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
        enkf_quit("OSISAF reader: could not convert file name \"%s\" to date", fname);
    basename[11] = 0;
    if (!str2int(&basename[35], &month))
        enkf_quit("OSISAF reader: could not convert file name \"%s\" to date", fname);
    basename[9] = 0;
    if (!str2int(&basename[35], &year))
        enkf_quit("OSISAF reader: could not convert file name \"%s\" to date", fname);
    snprintf(&tunits[tunits_len], MAXSTRLEN - tunits_len, " since %4d-%02d-%02d", year, month, day);

    ncw_close(ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    // load coast data from 'reference_grid.nc'
    filename="reference_grid.nc";
    ncw_open(filename, NC_NOWRITE, &ncid);
    // longitude
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "boundary_lon")) ? "boundary_lon" : "length", &lon_dimid);
    ncw_inq_dimlen(ncid, lon_dimid, &coast_len);
    ncw_inq_varid(ncid, "boundary_lon", &varid_lon);
    coast_lon = malloc(coast_len * sizeof(float));
    ncw_get_var_float(ncid, varid_lon, coast_lon);  
    // latitude
    ncw_inq_dimid(ncid, (ncw_dim_exists(ncid, "boundary_lat")) ? "boundary_lat" : "length", &lat_dimid);
    // ncw_inq_dimlen(ncid, lat_dimid, &coast_len);
    ncw_inq_varid(ncid, "boundary_lat", &varid_lat);
    coast_lat = malloc(coast_len * sizeof(float));
    ncw_get_var_float(ncid, varid_lat, coast_lat);
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
                o->instrument = st_add_ifabsent(obs->instruments, "osisaf", -1);
                o->id = obs->nobs;
                o->fid = fid;
                o->batch = 0;
                if (strcmp(meta->type,"SIC") == 0) {
                    o->value = (double) (sic[it][i][j]*sic_scale_factor*0.01);  // sic_scale_factor=0.01, sic[it][i][j]*sic_scale_factor is from 0 to 100.
                    o->estd   = 0.01+pow(0.5-fabs(0.5-o->value),2.0);
                }
                o->lon = lon[i][j];
                o->lat = lat[i][j];
                o->depth = 0.0;
                o->fk = (double) ksurf;
                o->status = grid_xy2fij_f(g, o->lon, o->lat, &o->fi, &o->fj);
                if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                    continue;
                // check distance between the observation the closest coast
                obs_distance2coast(g, coast_len, coast_lon, coast_lat, o->lon, o->lat, &o->status); 
                if (!obs->allobs && o->status == STATUS_SHALLOW)
                    continue;
                // if (!obs->allobs && o->status == STATUS_LAND)
                //     continue;
                o->model_depth = NAN;   // set in obs_add()
                o->time = time[tlen] * tunits_multiple + tunits_offset;
                o->aux = -1;

                obs->nobs++;
            }
        }
    }

    free(lon);
    free(lat);
    if (strcmp(meta->type,"SIC") == 0) 
        free(sic);
    
    free(time);
}

// /**
//  * @brief      filter out observations near coast (min_distance2coast_km)
//  * 
//  * @param g    grid object, containing coast coordinates.
//  * @param lon  longtitude of observation
//  * @param lat  latitude of observation
//  * @return int return status
//  */
// void obs_distance2coast(grid* g, size_t coast_len, float* coast_lon, float* coast_lat, float lon, float lat, short int *status){ 
//     int min_distance2coast_km =50;  //user defined distance between observation position and the closest coast. It is used to filter out observations too close to the coast.
//     double ll[2] = { lon, lat }; //observation coordinates
//     double xyz1[3];
//     double distance;
    
//     ll2xyz(ll, xyz1);
//     // distance filter
//     for (int j = 0; j < coast_len; ++j) {
//         double ll2[2] = { coast_lon[j], coast_lat[j] };
//         double xyz2[3];    
//         ll2xyz(ll2, xyz2);
//         // reject observation close to the coast, mark the obervation as STATUS_OUTSIDEGRID.
//         distance = sqrt((xyz1[0] - xyz2[0]) * (xyz1[0] - xyz2[0]) + (xyz1[1] -  xyz2[1]) * (xyz1[1] - xyz2[1]) + (xyz1[2] - xyz2[2]) * (xyz1[2] - xyz2[2]));
//         if (distance<=min_distance2coast_km){
//             *status = STATUS_SHALLOW; 
//             return;
//         }
//     }
//     // accept the observation.
//     return;
// }
