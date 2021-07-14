
/* Author: adpozuelo@gmail.com
 * Version: 1.0
 * Date: 03/2021
 */

#ifndef CONF_H
#define CONF_H

#include "mkl_vsl.h"

typedef float precision;
// typedef double precision;

#define INPUT_FILENAME "in.mc"
#define LOG_FILENAME "log.mc"
#define INPUT_XYZ_FILENAME "init.lammpstrj"
#define OUTPUT_XYZ_FILENAME "mc.lammpstrj"

#define NDIM 3
#define NTHREAD 64
#define MAX_LINE_SIZE 256

#define NSIM 10
#define OUPUT_TIME_FILENAME "nsim_times.mc"

#define CPU_TIME_INIT      \
  clock_t begin = clock(); \
  clock_t end = 0;

#define CPU_TIME_START begin = clock();

#define CPU_TIME_STOP \
  end = clock();      \
  cxf->time_spent += (double)(end - begin) / CLOCKS_PER_SEC;

#define CPU_TIME_STOP_MAIN \
  end = clock();      \
  cxf.time_spent += (double)(end - begin) / CLOCKS_PER_SEC;

typedef struct {
  VSLStreamStatePtr streamRNG;
  unsigned int nstep, nequil, natoms, ntrial, naccept;
  unsigned int *nspps;
  precision kt, acceptance, temp, density, sigma_o, volume, esr;
  char *mode, *units, *input_conf;
  char **atoms;
  unsigned short cuda_device, shift, nsp, nitmax, lammpstrj;
  unsigned short *ptype, *molecule;
  unsigned short **itp;
  precision *rdmax, *al, *bl, *bl2, *rc2, *r, *side, *esrrc;
  precision **rc;
  double time_spent;
} Configuration;

#endif