/* Author: adpozuelo@gmail.com
 * Version: 1.0
 * Date: 03/2021
 */

#ifndef IO_H
#define IO_H

#include "conf.h"

void read_input_file(Configuration *cxf);

void write_configuration(const Configuration *cxf, const unsigned int timestep);

void print_header(const Configuration *cxf);

void print_step(const Configuration *cxf, const unsigned int step);

#endif
