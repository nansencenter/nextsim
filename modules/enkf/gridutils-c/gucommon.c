/******************************************************************************
 *
 *  File:           gucommon.c
 *  
 *  Created         22/01/2002
 *  
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *  
 *  Purpose:        Some common stuff for grid utilities
 *  Revisions:      4 Jul 2013 PS
 *                    Using saved errno
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <assert.h>
#include <errno.h>
#include "version.h"
#include "gucommon.h"

#define BUFSIZE 10240

static void gu_quit_def(char* format, ...);

int gu_verbose = 0;
gu_quitfn gu_quit = gu_quit_def;

/**
 */
static void gu_quit_def(char* format, ...)
{
    va_list args;

    fflush(stdout);
    fprintf(stderr, "\n  error: gridutils: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");

    exit(1);
}

/**
 */
void gu_setquitfn(gu_quitfn quitfn)
{
    gu_quit = quitfn;
}

/**
 */
FILE* gu_fopen(const char* path, const char* mode)
{
    FILE* f = NULL;

    f = fopen(path, mode);
    if (f == NULL) {
        int errno_saved = errno;

        gu_quit("%s: could not open for \"%s\" : %s", path, mode, strerror(errno_saved));
    }

    return f;
}

/** Allocates ni x nj matrix of something and fills it with zeros. An element
 * (i,j) will be accessed as [j][i]. For deallocation use free().
 *
 * @param nj Dimension 2
 * @param ni Dimension 1
 * @param unitsize Size of one matrix element in bytes
 * @return Matrix
 */
void* gu_alloc2d(size_t nj, size_t ni, size_t unitsize)
{
    size_t size;
    void* p;
    void** pp;
    int i;

    if (ni <= 0 || nj <= 0)
        gu_quit("gu_alloc2d(): invalid size (nj = %d, ni = %d)", nj, ni);

    size = nj * sizeof(void*) + nj * ni * unitsize;
    if ((p = malloc(size)) == NULL) {
        int errno_saved = errno;

        gu_quit("alloc2d(): %s", strerror(errno_saved));
    }
    memset(p, 0, size);

    pp = p;
    p = &((size_t *) p)[nj];
    for (i = 0; i < nj; i++)
        pp[i] = &((char*) p)[i * ni * unitsize];

    return pp;
}

/** Destroys a matrix.
 * @param p Matrix
 */
void gu_free2d(void* p)
{
    free(p);
}

/**
 */
int** gu_readmask(char* fname, int nx, int ny)
{
    int** v = NULL;
    FILE* f = NULL;
    char buf[BUFSIZE];
    int count;
    int i, j;

    if (strcasecmp(fname, "stdin") == 0 || strcmp(fname, "-") == 0)
        f = stdin;
    else
        f = gu_fopen(fname, "r");

    v = gu_alloc2d(ny, nx, sizeof(int));

    for (j = 0, count = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            if (fgets(buf, BUFSIZE, f) == NULL)
                gu_quit("%s: could not read %d-th mask value (%d x %d values expected)", fname, j * nx + i + 1, nx, ny);
            buf[strlen(buf) - 1] = 0;
            if (strcmp(buf, "0") == 0)
                v[j][i] = 0;
            else if (strcmp(buf, "1") == 0) {
                v[j][i] = 1;
                count++;
            } else
                gu_quit("%s: could not interpret %d-th mask value = \"%s\" (expected \"0\" or \"1\"", fname, j * nx + i + 1, buf);
        }
    }

    if (f != stdin)
        fclose(f);

    if (gu_verbose) {
        int n = nx * ny;

        fprintf(stderr, "## mask: %d valid cells (%.1f%%), %d masked cells (%.1f%%)\n", count, 100.0 * count / n, n - count, 100.0 * (n - count) / n);
        fflush(stderr);
    }

    return v;
}
