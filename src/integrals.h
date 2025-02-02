/**************************************************************************
 *                                                                        *
 *   Author: Ivo Filot <ivo@ivofilot.nl>                                  *
 *                                                                        *
 *   P2000T-HF-H2 is free software:                                       *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   P2000T-HF-H2 is distributed in the hope that it will be useful,      *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#ifndef _INTEGRALS_H
#define _INTEGRALS_H

#include <math.h>
#include "memory.h"
#include "terminal.h"

#define PI 3.141527

float overlap_cgf(float *coeff1, float *coeff2, float *alpha1, float *alpha2, float *pos1, float *pos2, uint8_t sz1, uint8_t sz2);
float kinetic_cgf(float *coeff1, float *coeff2, float *alpha1, float *alpha2, float *pos1, float *pos2, uint8_t sz1, uint8_t sz2);
float kinetic(float alpha1, uint8_t l1, uint8_t m1, uint8_t n1, float ax, float ay, float az,
              float alpha2, uint8_t l2, uint8_t m2, uint8_t n2, float bx, float by, float bz);
float overlap_1d(uint8_t l1, uint8_t l2, float x1, float x2, float gamma);
float normalization_gto(float alpha, uint8_t l, uint8_t m, uint8_t n);
float overlap(float alpha1, uint8_t l1, uint8_t m1, uint8_t n1, float ax, float ay, float az,
              float alpha2, uint8_t l2, uint8_t m2, uint8_t n2, float bx, float by, float bz);
float binomial(int8_t a, int8_t b);
float binomial_prefactor(int8_t s, int8_t ia, int8_t ib, float xpa, float xpb);

#endif // __INTEGRALS_H