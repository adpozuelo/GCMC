/* Author: adpozuelo@gmail.com
 * Version: 1.0
 * Date: 03/2021
 */

#ifndef ENERGY_H
#define ENERGY_H

#include "conf.h"

precision lennard_jones(const precision r2, const short nit, const Configuration *cxf);

precision energy_cpu(const Configuration *cxf);

precision distance2(precision *r, const precision *side);

#endif
