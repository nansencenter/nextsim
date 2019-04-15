/******************************************************************************
 *
 * File:        obstypes.h       
 *
 * Created:     18/08/2014
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_OBSTYPES_H)

typedef struct {
    int id;
    char* name;
    int issurface;
    int nvar;
    char** varnames;
    char* alias;                /* variable name used in file names */
    char* offset_fname;
    char* offset_varname;
    char* mld_varname;
    double mld_threshold;
    char* hfunction;
    double allowed_min;
    double allowed_max;
    int isasync;
    double async_tstep;
    int async_centred;
    int nlocrad;
    double* locrad;
    double* locweight;
    double rfactor;
    int nlobsmax;
    double estdmin;

    int vid;
    int gridid;
    int sob_stride;

    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;

    int nobs;
    int ngood;
    int noutside_grid;
    int noutside_obsdomain;
    int noutside_obswindow;
    int nland;
    int nshallow;
    int nbadbatch;
    int nrange;
    int nsubgrid;
    int nthinned;
    int nmodified;

    /*
     * allowed time range of observations 
     */
    double windowmin;
    double windowmax;

    /*
     * actual time range of observations 
     */
    double date_min;
    double date_max;

    /*
     * domains observations of this type are visible from
     */
    int ndomains;
    char** domainnames;
} obstype;

void obstypes_read(enkfprm* prm, char fname[], int* n, obstype** types);
void obstypes_describeprm(void);
void obstypes_set(int n, obstype* types, model* m);
void obstypes_destroy(int n, obstype* types);

int obstype_getid(int n, obstype types[], char* name, int hastosucceed);
double obstype_calclcoeff(obstype* type, double dist);
double obstype_getmaxlocrad(obstype* type);

#define _OBSTYPES_H
#endif
