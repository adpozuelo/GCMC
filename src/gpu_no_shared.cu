/* Author: adpozuelo@gmail.com
 * Version: 1.3
 * Date: 07/2021
 */

#include <stdio.h>
#include <stdlib.h>

#include "gpu.h"

extern "C" {
#include "conf.h"
#include "energy.h"
#include "mkl_vsl.h"
}

__device__ precision __distance2__(precision *r, const precision *side) {
  for (int i = 0; i < NDIM; ++i) {
    if (r[i] > 0.5) r[i] -= 1;
    if (r[i] < -0.5) r[i] += 1;
  }

  precision rd2 = 0.0;
  for (int i = 0; i < NDIM; ++i) {
    r[i] *= side[i];
    r[i] *= r[i];
    rd2 += r[i];
  }

  return rd2;
}

__device__ precision __lennard_jones__(const precision r2,
                                       const unsigned short nit,
                                       const precision *al,
                                       const precision *bl2) {
  precision r6 = (bl2[nit] / r2) * (bl2[nit] / r2) * (bl2[nit] / r2);
  return 4 * al[nit] * r6 * (r6 - 1.0);
}

__global__ void __binary_reduction__(unsigned int *natoms_nsp,
                                     precision *g_idata, precision *g_odata) {
  __shared__ precision sdata[NTHREAD];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < natoms_nsp[0])
    sdata[tid] = g_idata[i];
  else
    sdata[tid] = 0.0;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) sdata[tid] += sdata[tid + s];
    __syncthreads();
  }

  if (tid == 0) atomicAdd(g_odata, sdata[0]);
}

__global__ void energy_gpu(precision *eng, unsigned int *natoms_nsp,
                           unsigned short *itp, unsigned short *ptype,
                           precision *r, precision *side, precision *rc2,
                           precision *al, precision *bl2, precision *esrrc) {
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int natoms = natoms_nsp[0];
  unsigned int nsp = natoms_nsp[1];
  precision energ = 0.0;

  if (i < natoms) {
    precision rd2;
    precision rdd[NDIM];
    unsigned short nit;

    unsigned int imol = i * NDIM;
    precision xi = r[imol];
    precision yi = r[imol + 1];
    precision zi = r[imol + 2];

    for (unsigned int j = 0; j < natoms; j++) {
      if (i != j) {
        unsigned int jmol = j * NDIM;
        nit = itp[ptype[i] * nsp + ptype[j]];

        rdd[0] = xi - r[jmol];
        rdd[1] = yi - r[jmol + 1];
        rdd[2] = zi - r[jmol + 2];

        rd2 = __distance2__(rdd, side);
        if (rd2 < rc2[nit])
          energ += __lennard_jones__(rd2, nit, al, bl2) - esrrc[nit];
      }
    }
  }
  eng[i] = energ;
}

__global__ void delta_energy_gpu(unsigned int ntest, precision *eng0,
                                 precision *eng1, precision *r_test,
                                 unsigned int *natoms_nsp, unsigned short *itp,
                                 unsigned short *ptype, precision *r,
                                 precision *side, precision *rc2, precision *al,
                                 precision *bl2, precision *esrrc) {
  unsigned int j = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int natoms = natoms_nsp[0];
  unsigned int nsp = natoms_nsp[1];

  precision energ0 = 0.0;
  precision energ1 = 0.0;

  if (j < natoms) {
    precision rdd[NDIM], rddn[NDIM];
    precision rd2;
    unsigned short nit;
    if (j != ntest) {
      unsigned int jmol = j * NDIM;
      unsigned int test_mol = ntest * NDIM;
      for (int k = 0; k < NDIM; ++k) {
        rdd[k] = r[jmol + k] - r[test_mol + k];
        rddn[k] = r[jmol + k] - r_test[k];
      }

      nit = itp[ptype[ntest] * nsp + ptype[j]];

      // before movement
      rd2 = __distance2__(rdd, side);
      if (rd2 < rc2[nit])
        energ0 = __lennard_jones__(rd2, nit, al, bl2) - esrrc[nit];

      // after movement
      rd2 = __distance2__(rddn, side);
      if (rd2 < rc2[nit])
        energ1 = __lennard_jones__(rd2, nit, al, bl2) - esrrc[nit];
    }
  }
  eng0[j] = energ0;
  eng1[j] = energ1;
}

extern "C" void gpu(Configuration *cxf, const unsigned short mode) {
  CPU_TIME_INIT
  GPU_TIME_INIT

  static unsigned int *natoms_nsp_dev;
  static unsigned short *itp_dev, *ptype_dev;
  static precision *r_dev, *side_dev, *rc2_dev, *al_dev, *bl2_dev, *esrrc_dev;

  if (mode == 0) {  // Initialize GPU memory
    unsigned short nsp2 = cxf->nsp * cxf->nsp;
    unsigned short *itp_serialized =
        (unsigned short *)malloc(nsp2 * sizeof(unsigned short));
    for (int i = 0; i < cxf->nsp; ++i)
      for (int j = 0; j < cxf->nsp; ++j)
        itp_serialized[i * cxf->nsp + j] = cxf->itp[i][j];
    unsigned int natoms_nsp[2] = {cxf->natoms, cxf->nsp};
    cudaSetDevice(cxf->cuda_device);
    cudaMalloc((void **)&natoms_nsp_dev, 2 * sizeof(unsigned int));
    cudaMalloc((void **)&itp_dev, nsp2 * sizeof(unsigned short));
    cudaMalloc((void **)&ptype_dev, cxf->natoms * sizeof(unsigned short));
    cudaMalloc((void **)&r_dev, cxf->natoms * NDIM * sizeof(precision));
    cudaMalloc((void **)&side_dev, NDIM * sizeof(precision));
    cudaMalloc((void **)&rc2_dev, cxf->nitmax * sizeof(precision));
    cudaMalloc((void **)&al_dev, cxf->nitmax * sizeof(precision));
    cudaMalloc((void **)&bl2_dev, cxf->nitmax * sizeof(precision));
    cudaMalloc((void **)&esrrc_dev, cxf->nitmax * sizeof(precision));

    cudaMemcpy(natoms_nsp_dev, natoms_nsp, 2 * sizeof(unsigned int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(itp_dev, itp_serialized, nsp2 * sizeof(unsigned short),
               cudaMemcpyHostToDevice);
    cudaMemcpy(ptype_dev, cxf->ptype, cxf->natoms * sizeof(unsigned short),
               cudaMemcpyHostToDevice);
    cudaMemcpy(r_dev, cxf->r, cxf->natoms * NDIM * sizeof(precision),
               cudaMemcpyHostToDevice);
    cudaMemcpy(side_dev, cxf->side, NDIM * sizeof(precision),
               cudaMemcpyHostToDevice);
    cudaMemcpy(rc2_dev, cxf->rc2, cxf->nitmax * sizeof(precision),
               cudaMemcpyHostToDevice);
    cudaMemcpy(al_dev, cxf->al, cxf->nitmax * sizeof(precision),
               cudaMemcpyHostToDevice);
    cudaMemcpy(bl2_dev, cxf->bl2, cxf->nitmax * sizeof(precision),
               cudaMemcpyHostToDevice);
    cudaMemcpy(esrrc_dev, cxf->esrrc, cxf->nitmax * sizeof(precision),
               cudaMemcpyHostToDevice);
    CPU_TIME_STOP
  } else if (mode == 1) {  // energy_gpu
    precision *e_by_thread_dev, *total_esr_dev, total_esr;
    cudaMalloc((void **)&e_by_thread_dev, cxf->natoms * sizeof(precision));
    cudaMalloc((void **)&total_esr_dev, sizeof(precision));

    unsigned int nblock = cxf->natoms / NTHREAD;
    if (cxf->natoms % NTHREAD != 0) ++nblock;

    cudaMemset(total_esr_dev, 0, sizeof(precision));
    CPU_TIME_STOP

    GPU_TIME_START
    energy_gpu<<<nblock, NTHREAD>>>(e_by_thread_dev, natoms_nsp_dev, itp_dev,
                                    ptype_dev, r_dev, side_dev, rc2_dev, al_dev,
                                    bl2_dev, esrrc_dev);

    __binary_reduction__<<<nblock, NTHREAD>>>(natoms_nsp_dev, e_by_thread_dev,
                                              total_esr_dev);
    GPU_TIME_STOP

    CPU_TIME_START
    cudaMemcpy(&total_esr, total_esr_dev, sizeof(precision),
               cudaMemcpyDeviceToHost);

    cxf->esr = total_esr / 2;
    cudaFree(e_by_thread_dev);
    cudaFree(total_esr_dev);
    CPU_TIME_STOP
  } else if (mode == 2) {  // move_atoms_gpu
    unsigned int ntest;
    const int harvest_size = NDIM + 1;
    precision deltae, e_before, e_after;
    double *harvest = (double *)malloc(harvest_size * sizeof(double));

    precision *r_test = (precision *)malloc(NDIM * sizeof(precision));
    precision *r_test_dev;
    cudaMalloc((void **)&r_test_dev, NDIM * sizeof(precision));

    precision *e_before_by_thread_dev, *e_after_by_thread_dev, *e_before_dev,
        *e_after_dev;
    cudaMalloc((void **)&e_before_by_thread_dev,
               cxf->natoms * sizeof(precision));
    cudaMalloc((void **)&e_after_by_thread_dev,
               cxf->natoms * sizeof(precision));
    cudaMalloc((void **)&e_before_dev, sizeof(precision));
    cudaMalloc((void **)&e_after_dev, sizeof(precision));

    unsigned int nblock = (cxf->natoms + (NTHREAD - 1)) / NTHREAD;

    for (int i = 0; i < cxf->natoms; ++i) {
      cxf->ntrial++;
      vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, cxf->streamRNG, harvest_size,
                   harvest, 0, 1);
      ntest = (unsigned int)cxf->natoms * harvest[NDIM];

      for (int j = 0; j < NDIM; ++j) {
        r_test[j] = cxf->r[ntest * NDIM + j] +
                    cxf->rdmax[j] * (2 * harvest[j] - 1) / cxf->side[j];
        if (r_test[j] < 0) r_test[j] += 1;
        if (r_test[j] > 1) r_test[j] -= 1;
      }

      cudaMemcpy(r_test_dev, r_test, NDIM * sizeof(precision),
                 cudaMemcpyHostToDevice);
      cudaMemset(e_before_dev, 0, sizeof(precision));
      cudaMemset(e_after_dev, 0, sizeof(precision));
      CPU_TIME_STOP

      GPU_TIME_START
      delta_energy_gpu<<<nblock, NTHREAD>>>(
          ntest, e_before_by_thread_dev, e_after_by_thread_dev, r_test_dev,
          natoms_nsp_dev, itp_dev, ptype_dev, r_dev, side_dev, rc2_dev, al_dev,
          bl2_dev, esrrc_dev);

      __binary_reduction__<<<nblock, NTHREAD>>>(
          natoms_nsp_dev, e_before_by_thread_dev, e_before_dev);

      __binary_reduction__<<<nblock, NTHREAD>>>(
          natoms_nsp_dev, e_after_by_thread_dev, e_after_dev);
      GPU_TIME_STOP

      CPU_TIME_START
      cudaMemcpy(&e_before, e_before_dev, sizeof(precision),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&e_after, e_after_dev, sizeof(precision),
                 cudaMemcpyDeviceToHost);

      deltae = e_after - e_before;

      if (deltae < 0.0) {
        for (int k = 0; k < NDIM; ++k) cxf->r[ntest * NDIM + k] = r_test[k];

        cudaMemcpy(r_dev + ntest * NDIM, r_test, NDIM * sizeof(precision),
                   cudaMemcpyHostToDevice);

        cxf->esr += deltae;
        cxf->naccept++;
      } else {
        double xi[1];
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, cxf->streamRNG, 1, xi, 0, 1);
        if (exp(-deltae) > xi[0]) {
          for (int k = 0; k < NDIM; ++k) cxf->r[ntest * NDIM + k] = r_test[k];

          cudaMemcpy(r_dev + ntest * NDIM, r_test, NDIM * sizeof(precision),
                     cudaMemcpyHostToDevice);

          cxf->esr += deltae;
          cxf->naccept++;
        }
      }
    }

    cudaFree(e_before_by_thread_dev);
    cudaFree(e_after_by_thread_dev);
    cudaFree(r_test_dev);
    cudaFree(e_before_dev);
    cudaFree(e_after_dev);
    free(harvest);
    free(r_test);
    CPU_TIME_STOP
  } else if (mode == 3) {  // Release GPU memory
    cudaFree(natoms_nsp_dev);
    cudaFree(itp_dev);
    cudaFree(ptype_dev);
    cudaFree(r_dev);
    cudaFree(side_dev);
    cudaFree(rc2_dev);
    cudaFree(al_dev);
    cudaFree(bl2_dev);
    cudaFree(esrrc_dev);
    CPU_TIME_STOP
  } else {
    fputs("ERROR: Incorrect GPU code!\n", stderr);
    exit(1);
  }
  GPU_TIME_DESTROY
}
