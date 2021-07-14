/* Author: adpozuelo@gmail.com
 * Version: 1.0
 * Date: 03/2021
 */

#ifndef GPU_H
#define GPU_H

#include "conf.h"

#define GPU_TIME_INIT      \
  cudaEvent_t start, stop; \
  cudaEventCreate(&start); \
  cudaEventCreate(&stop);  \
  float elapsed_time;

#define GPU_TIME_DESTROY   \
  cudaEventDestroy(start); \
  cudaEventDestroy(stop);

#define GPU_TIME_START cudaEventRecord(start);

#define GPU_TIME_STOP                               \
  cudaEventRecord(stop);                            \
  cudaEventSynchronize(stop);                       \
  elapsed_time = 0.0;                               \
  cudaEventElapsedTime(&elapsed_time, start, stop); \
  cxf->time_spent += elapsed_time / 1000.0;

extern "C" void gpu(Configuration *cxf, const unsigned short mode);

#endif