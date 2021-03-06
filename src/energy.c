/* Author: adpozuelo@gmail.com
 * Version: 1.0
 * Date: 03/2021
 */

#include <stdlib.h>
#include "energy.h"

precision distance2(precision *r, const precision *side)
{
    for (int i = 0; i < NDIM; ++i)
    {
        if (r[i] > 0.5)
            r[i] -= 1;
        if (r[i] < -0.5)
            r[i] += 1;
    }

    precision rd2 = 0.0;
    for (int i = 0; i < NDIM; ++i)
    {
        r[i] *= side[i];
        r[i] *= r[i];
        rd2 += r[i];
    }

    return rd2;
}

precision lennard_jones(const precision r2, const short nit, const Configuration *cxf)
{
    precision r6 = (cxf->bl2[nit] / r2) * (cxf->bl2[nit] / r2) * (cxf->bl2[nit] / r2);
    return 4 * cxf->al[nit] * r6 * (r6 - 1.0);
}

precision energy_cpu(const Configuration *cxf)
{
    precision rd2, energy = 0.0;
    precision *rdd = (precision *)malloc(NDIM * sizeof(precision));
    unsigned short nit;

    for (unsigned int i = 0; i < cxf->natoms; ++i)
    {
        for (unsigned int j = i + 1; j < cxf->natoms; ++j)
        {
            nit = cxf->itp[cxf->ptype[i]][cxf->ptype[j]];
            for (int k = 0; k < NDIM; ++k)
                rdd[k] = cxf->r[i * NDIM + k] - cxf->r[j * NDIM + k];
            rd2 = distance2(rdd, cxf->side);
            if (rd2 < cxf->rc2[nit])
                energy += lennard_jones(rd2, nit, cxf) - cxf->esrrc[nit];
        }
    }
    
    free(rdd);
    return energy;
}