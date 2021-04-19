/******************************************************************************
 *
 * File:           grikmap.c
 *  
 * Created:        16 March 2016
 *  
 * Author:         Pavel Sakov
 *                 BoM
 *
 * Purpose:        Mapping of curvilinear grids based on kd-tree.
 *
 * Revisions:      
 *
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include "kdtree.h"
#include "poly.h"
#include "gridkmap.h"
#include "gucommon.h"

#define EPS 1.0e-8
#define EPS_ZERO 1.0e-5

#define _REARTH 6371.0
#define _DEG2RAD (M_PI / 180.0)

struct gridkmap {
    int nce1;                   /* number of cells in e1 direction */
    int nce2;                   /* number of cells in e2 direction */
    double** gx;                /* reference to array of X coords
                                 * [nce2+1][nce1+1] */
    double** gy;                /* reference to array of Y coords
                                 * [nce2+1][nce1+1] */
    double** gX; 
    double** gY;
    double** gZ;
    kdtree* tree;               /* kd tree with grid nodes */
};

/** Builds a grid map structure to facilitate conversion from coordinate
 * to index space.
 *
 * @param gx array of X coordinates [nce2 + 1][nce1 + 1]
 * @param gy array of Y coordinates [nce2 + 1][nce1 + 1]
 * @param nce1 number of cells in e1 direction
 * @param nce2 number of cells in e2 direction
 * @return a map tree to be used by xy2ij
 */
gridkmap* gridkmap_build(int nce1, int nce2, double** gx, double** gy)
{
    gridkmap* gm = malloc(sizeof(gridkmap));
    double* data[3];
    double xyz[3];
    double ll[2];
    int i, j;
    gm->gX = gu_alloc2d(nce2+1, nce1+1, sizeof(double));
    gm->gY = gu_alloc2d(nce2+1, nce1+1, sizeof(double));
    gm->gZ = gu_alloc2d(nce2+1, nce1+1, sizeof(double));

    gm->nce1 = nce1;
    gm->nce2 = nce2;
    gm->gx = gx;
    gm->gy = gy;

    gm->tree = kd_create(3);
    // data[0] = gx[0];
    // data[1] = gy[0];
    for (j=0; j<nce2+1; ++j){
        for (i=0; i<nce1+1; ++i){
            ll[0] = gx[j][i];
            ll[1] = gy[j][i];
            _ll2xyz(ll, xyz);
            gm->gX[j][i] = xyz[0]; gm->gY[j][i] = xyz[1]; gm->gZ[j][i] = xyz[2]; 
        }
    }
    data[0] = gm->gX[0];
    data[1] = gm->gY[0];
    data[2] = gm->gZ[0];
    kd_insertnodes(gm->tree, (nce1 + 1) * (nce2 + 1), data, 1 /* shuffle */ );
    // kd_insertnodes(gm->tree, (nce1 + 1) * (nce2 + 1), data, 1 /* shuffle */ );

    return gm;
}

/**
 */
void gridkmap_destroy(gridkmap* gm)
{
    kd_destroy(gm->tree);
    free(gm);
}

/**
 */
int gridkmap_xy2ij(gridkmap* gm, double x, double y, int* iout, int* jout)
{
    double pos[3];
    double ll[2];
    size_t nearest;
    size_t id;
    poly* p;
    int i, j, i1, i2, j1, j2;
    int success = 0;


    ll[0] = x;
    ll[1] = y;

    _ll2xyz(ll, pos);

    double* minmax = kd_getminmax(gm->tree);
    if (pos[0] < minmax[0] || pos[1] < minmax[1] || pos[0] > minmax[2] || pos[1] > minmax[3])
        return success;

    nearest = kd_findnearestnode(gm->tree, pos);
    id = kd_getnodeorigid(gm->tree, nearest);
    p = poly_create();

    j = id / (gm->nce1 + 1);
    i = id % (gm->nce1 + 1);

    i1 = (i > 0) ? i - 1 : i;
    i2 = (i < gm->nce1) ? i + 1 : i;
    j1 = (j > 0) ? j - 1 : j;
    j2 = (j < gm->nce2) ? j + 1 : j;

    /*
     * TODO?: looks a bit heavyweight
     */
    for (j = j1; j <= j2 - 1; ++j)
        for (i = i1; i <= i2 - 1; ++i) {
            if (!isfinite(gm->gx[j][i]) || !isfinite(gm->gx[j][i + 1]) || !isfinite(gm->gx[j + 1][i + 1]) || !isfinite(gm->gx[j + 1][i]))
                continue;
            poly_addpoint(p, gm->gx[j][i], gm->gy[j][i]);
            poly_addpoint(p, gm->gx[j][i + 1], gm->gy[j][i + 1]);
            poly_addpoint(p, gm->gx[j + 1][i + 1], gm->gy[j + 1][i + 1]);
            poly_addpoint(p, gm->gx[j + 1][i], gm->gy[j + 1][i]);
            poly_addpoint(p, gm->gx[j][i], gm->gy[j][i]);

            if (poly_containspoint(p, x, y)) {
                success = 1;
                *iout = i;
                *jout = j;
                goto finish;
            }
            poly_clear(p);
        }

  finish:
    poly_destroy(p);

    return success;
}

/**
 */
int gridkmap_getnce1(gridkmap* gm)
{
    return gm->nce1;
}

/**
 */
int gridkmap_getnce2(gridkmap* gm)
{
    return gm->nce2;
}

/**
 */
double** gridkmap_getxnodes(gridkmap* gm)
{
    return gm->gx;
}

/**
 */
double** gridkmap_getynodes(gridkmap* gm)
{
    return gm->gy;
}

/** Converts from geographic to cartesian coordinates.
 * @param in Input: {lon, lat}
 * @param out Output: {x, y, z}
 */
void _ll2xyz(double in[2], double* out)
{
    double lon = in[0] * _DEG2RAD;
    double lat = in[1] * _DEG2RAD;
    double coslat = cos(lat);

    out[0] = _REARTH * sin(lon) * coslat;
    out[1] = _REARTH * cos(lon) * coslat;
    out[2] = _REARTH * sin(lat);
}
