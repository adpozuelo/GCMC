/* Author: adpozuelo@gmail.com
 * Version: 1.0
 * Date: 03/2021
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "xyz.h"
#include "mkl_vsl.h"

void init_fcc(Configuration *cxf)
{
    precision ratom[4][3];
    precision boxl = pow((precision)cxf->natoms / cxf->density, 0.33333333333) / cxf->sigma_o;

    cxf->side = (precision *)malloc(NDIM * sizeof(precision));
    cxf->volume = 1.0;
    for (int i = 0; i < NDIM; ++i)
    {
        cxf->side[i] = boxl;
        cxf->volume *= cxf->side[i];
    }

    cxf->r = (precision *)malloc(cxf->natoms * NDIM * sizeof(precision));

    unsigned int mside = (unsigned int)ceil(pow((double)cxf->natoms / 4.0, 0.333333));
    precision boxfc2 = boxl / mside;

    for (int i = 0; i < NDIM; ++i)
    {
        ratom[0][i] = 0;
        ratom[1][i] = boxfc2 / 2;
        ratom[2][i] = boxfc2 / 2;
        ratom[3][i] = boxfc2 / 2;
    }
    ratom[1][0] = 0;
    ratom[2][1] = 0;
    ratom[3][2] = 0;

    unsigned int m = 0;
    for (unsigned int i = 0; i < mside; ++i)
    {
        for (unsigned int j = 0; j < mside; ++j)
        {
            for (unsigned int k = 0; k < mside; ++k)
            {
                precision despx = i * boxfc2;
                precision despy = j * boxfc2;
                precision despz = k * boxfc2;
                for (int l = 0; l < 4; ++l)
                {
                    if (m < 3 * cxf->natoms)
                    {
                        cxf->r[m] = ratom[l][0] + despx - boxl / 2;
                        cxf->r[m + 1] = ratom[l][1] + despy - boxl / 2;
                        cxf->r[m + 2] = ratom[l][2] + despz - boxl / 2;
                        m += 3;
                    }
                }
            }
        }
    }

    unsigned short shuffle_times = 2;
    precision xtmp, ytmp, ztmp;
    int *random1 = (int *)malloc(cxf->natoms * sizeof(int));
    int *random2 = (int *)malloc(cxf->natoms * sizeof(int));
    for (int p = 0; p < shuffle_times; ++p)
    {
        viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, cxf->streamRNG, cxf->natoms, random1, 0, cxf->natoms);
        viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, cxf->streamRNG, cxf->natoms, random2, 0, cxf->natoms);
        for (unsigned int i = 0; i < cxf->natoms; ++i)
        {
            xtmp = cxf->r[random1[i] * NDIM];
            ytmp = cxf->r[random1[i] * NDIM + 1];
            ztmp = cxf->r[random1[i] * NDIM + 2];
            cxf->r[random1[i] * NDIM] = cxf->r[random2[i] * NDIM];
            cxf->r[random1[i] * NDIM + 1] = cxf->r[random2[i] * NDIM + 1];
            cxf->r[random1[i] * NDIM + 2] = cxf->r[random2[i] * NDIM + 2];
            cxf->r[random2[i] * NDIM] = xtmp;
            cxf->r[random2[i] * NDIM + 1] = ytmp;
            cxf->r[random2[i] * NDIM + 2] = ztmp;
        }
    }
    free(random1);
    free(random2);

    for (int i = 1; i < cxf->nsp; ++i)
        for (int j = 0; j < i; ++j)
            cxf->nspps[i] += cxf->nspps[j];

    cxf->ptype = (unsigned short *)malloc(cxf->natoms * sizeof(unsigned short));
    cxf->molecule = (unsigned short *)calloc(cxf->natoms, sizeof(unsigned short));
    int x = 0;
    for (int i = 0; i < cxf->nsp; ++i)
    {
        for (; x < cxf->nspps[i]; ++x)
        {
            cxf->ptype[x] = i;
        }
    }

    for (unsigned int i = 0; i < cxf->natoms; ++i)
        for (int j = 0; j < NDIM; ++j)
        {
            cxf->r[i * NDIM + j] += cxf->side[j] / 2.0;
            cxf->r[i * NDIM + j] /= cxf->side[j];
        }
}

void read_lammpstrj(Configuration *cxf)
{
    char line[MAX_LINE_SIZE];
    FILE *file = fopen(INPUT_XYZ_FILENAME, "r");

    fgets(line, MAX_LINE_SIZE, file);
    fgets(line, MAX_LINE_SIZE, file);
    unsigned short timestep;
    unsigned int natoms;
    sscanf(line, "%hu", &timestep);
    if (timestep != 0)
    {
        printf("ERROR! no first timestep!\n");
        exit(1);
    }

    fgets(line, MAX_LINE_SIZE, file);
    fgets(line, MAX_LINE_SIZE, file);
    sscanf(line, "%u", &natoms);
    if (natoms != cxf->natoms)
    {
        printf("ERROR! natoms not equals!\n");
        exit(1);
    }

    fgets(line, MAX_LINE_SIZE, file);
    double trash_p, side;
    cxf->side = (precision *)malloc(NDIM * sizeof(precision));
    cxf->volume = 1.0;
    for (int i = 0; i < NDIM; ++i)
    {
        fgets(line, MAX_LINE_SIZE, file);
        sscanf(line, "%le %le", &trash_p, &side);
        cxf->side[i] = (precision)side;
        cxf->side[i] /= cxf->sigma_o;
        cxf->volume *= cxf->side[i];
    }

    unsigned int trash_ui;
    double r1, r2, r3;
    cxf->r = (precision *)malloc(cxf->natoms * NDIM * sizeof(precision));
    cxf->ptype = (unsigned short *)malloc(cxf->natoms * sizeof(unsigned short));
    cxf->molecule = (unsigned short *)malloc(cxf->natoms * sizeof(unsigned short));
    fgets(line, MAX_LINE_SIZE, file);
    for (unsigned int i = 0; i < cxf->natoms; ++i)
    {
        fgets(line, MAX_LINE_SIZE, file);
        sscanf(line, "%u %hu %hu %lf %lf %lf", &trash_ui, &cxf->ptype[i], &cxf->molecule[i], &r1, &r2, &r3);
        cxf->r[i * NDIM] = (precision)r1;
        cxf->r[i * NDIM + 1] = (precision)r2;
        cxf->r[i * NDIM + 2] = (precision)r3;
    }

    fclose(file);

    for (unsigned int i = 0; i < cxf->natoms; ++i)
    {
        cxf->ptype[i] -= 1;
        for (int j = 0; j < NDIM; ++j)
        {
            cxf->r[i * NDIM + j] /= cxf->side[j];
            cxf->r[i * NDIM + j] /= cxf->sigma_o;
        }
    }

    for (int i = 1; i < cxf->nsp; ++i)
        for (int j = 0; j < i; ++j)
            cxf->nspps[i] += cxf->nspps[j];
}
