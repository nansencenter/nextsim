/******************************************************************************
 *
 * File:        allreaders.c        
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stringtable.h"
#include "kdtree.h"
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

typedef struct {
    char* product;
    char* reader;
    obsread_fn readfn;
} obsreader_entry;

obsreader_entry allreaders[] = {
    {"NAVO", "standard", reader_navo_standard},
    {"WINDSAT", "standard", reader_windsat_standard},
    {"CARS", "standard", reader_cars_standard},
    {"MMT", "standard", reader_mmt_standard},
    {"AMSR2", "standard", reader_amsr2_standard},
    {"AMSRE", "standard", reader_amsre_standard},
    {"HIMAWARI8", "standard", reader_h8_standard},
    {"CMEMS", "standard", reader_cmems_standard},
    {"ANY", "scattered", reader_scattered},
    {"ANY", "xy_gridded", reader_xy_gridded},
    {"ANY", "xyz_gridded", reader_xyz_gridded},
    {"ANY", "xyh_gridded", reader_xyh_gridded},
    {"ANY", "z", reader_z},
    {"CS2SMOS", "standard", reader_cs2smos_standard},
};

/**
 */
static void describe_readers(void)
{
    int nreaders = sizeof(allreaders) / sizeof(obsreader_entry);
    int i;

    enkf_printf("  Available readers:\n");
    enkf_printf("    product   reader\n");
    enkf_printf("    ----------------\n");
    for (i = 0; i < nreaders; ++i)
        enkf_printf("    %-8s  %s\n", allreaders[i].product, allreaders[i].reader);
}

/**
 */
obsread_fn get_obsreadfn(obsmeta* m)
{
    int nreaders = sizeof(allreaders) / sizeof(obsreader_entry);
    int i;

    for (i = 0; i < nreaders; ++i)
        if ((strcmp(allreaders[i].product, "ANY") == 0 || strcasecmp(m->product, allreaders[i].product) == 0) && strcmp(m->reader, allreaders[i].reader) == 0)
            return allreaders[i].readfn;

    enkf_printf("\n\n  ERROR: no observation reader \"%s\" for product \"%s\"\n\n", m->reader, m->product);
    describe_readers();
    enkf_quit("bailing out", m->product);
    return NULL;
}
