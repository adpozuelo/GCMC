/* Author: adpozuelo@gmail.com
 * Version: 1.1
 * Date: 05/2021
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "moves.h"
#include "energy.h"

void move_atoms_cpu(Configuration *cxf)
{
    unsigned int ntest;
    unsigned short nit;
    const int harvest_size = NDIM + 1;
    precision distance2_before, distance2_after, delta_energy, energy_before, energy_after;
    double *harvest = (double *)malloc(harvest_size * sizeof(double));
    precision *r_test = (precision *)malloc(NDIM * sizeof(precision));
    precision *distance_before = (precision *)malloc(NDIM * sizeof(precision));
    precision *distance_after = (precision *)malloc(NDIM * sizeof(precision));

    for (unsigned int i = 0; i < cxf->natoms; ++i)
    {
        cxf->ntrial++;
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, cxf->streamRNG, harvest_size, harvest, 0, 1);
        ntest = (unsigned int)cxf->natoms * harvest[NDIM];

        for (int j = 0; j < NDIM; ++j)
        {
            r_test[j] = cxf->r[ntest * NDIM + j] + cxf->rdmax[j] * (2 * harvest[j] - 1) / cxf->side[j];
            if (r_test[j] < 0)
                r_test[j] += 1;
            if (r_test[j] > 1)
                r_test[j] -= 1;
        }

        energy_before = energy_after = 0.0;
        for (unsigned int j = 0; j < cxf->natoms; ++j)
        {
            if (j == ntest)
                continue;

            nit = cxf->itp[cxf->ptype[ntest]][cxf->ptype[j]];

            for (int k = 0; k < NDIM; ++k)
            {
                distance_before[k] = cxf->r[ntest * NDIM + k] - cxf->r[j * NDIM + k];
                distance_after[k] =  r_test[k] - cxf->r[j * NDIM + k];
            }

            distance2_before = distance2(distance_before, cxf->side);
            if (distance2_before < cxf->rc2[nit])
                energy_before += lennard_jones(distance2_before, nit, cxf) - cxf->esrrc[nit];

            distance2_after = distance2(distance_after, cxf->side);
            if (distance2_after < cxf->rc2[nit])
                energy_after += lennard_jones(distance2_after, nit, cxf) - cxf->esrrc[nit];
        }

        delta_energy = energy_after - energy_before;
        if (delta_energy < 0.0)
        {
            for (int k = 0; k < NDIM; ++k)
                cxf->r[ntest * NDIM + k] = r_test[k];

            cxf->esr += delta_energy;
            cxf->naccept++;
        }
        else
        {
            double xi[1];
            vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, cxf->streamRNG, 1, xi, 0, 1);

            if (exp(-delta_energy) > xi[0])
            {
                for (int k = 0; k < NDIM; ++k)
                    cxf->r[ntest * NDIM + k] = r_test[k];

                cxf->esr += delta_energy;
                cxf->naccept++;
            }
        }
    }
    
    free(harvest);
    free(r_test);
    free(distance_before);
    free(distance_after);
}