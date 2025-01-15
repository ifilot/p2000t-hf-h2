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

#include "integrals.h"

// pre-define factorials
float factorial[] = {1,1,2,6,24,120,720,5040,40320,362880,3628800};
float double_factorial[] = {1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840};

float overlap_cgf(float *coeff1, float *coeff2, float *alpha1, float *alpha2, 
                  float *pos1, float *pos2, uint8_t sz1, uint8_t sz2) {
    float s = 0.0;
    for(uint8_t i=0; i<sz1; i++) {
        for(uint8_t j=0; j<sz2; j++) {
            s += coeff1[i] * coeff2[j] * overlap(alpha1[i], 0, 0, 0, pos1[0], pos1[1], pos1[2], 
                                                 alpha2[j], 0, 0, 0, pos2[0], pos2[1], pos2[2]);
        }
    }

    return s;
}

float overlap(float alpha1, uint8_t l1, uint8_t m1, uint8_t n1, float ax, float ay, float az,
              float alpha2, uint8_t l2, uint8_t m2, uint8_t n2, float bx, float by, float bz) {

    float sum = 0.0;

    float rab2 = sqrt((ax-bx)*(ax-bx) + (ay-by)*(ay-by) + (az-bz)*(az-bz));
    float gamma = alpha1 + alpha2;
    float px = alpha1 * ax + alpha2 * bx / (alpha1 + alpha2);
    float py = alpha1 * ay + alpha2 * by / (alpha1 + alpha2);
    float pz = alpha1 * az + alpha2 * bz / (alpha1 + alpha2);

    float pre = pow(PI / gamma, 1.5) * exp(-alpha1 * alpha2 * rab2 / gamma);
    float wx = overlap_1d(l1, l2, px-ax, px-bx, gamma);
    float wy = overlap_1d(m1, m2, py-ay, py-by, gamma);
    float wz = overlap_1d(n1, n2, pz-az, pz-bz, gamma);

    // sprintf(termbuffer, "%6.4f %6.4f %6.4f %6.4f", pre, wx, wy, wz);
    // terminal_printtermbuffer();

    return pre * wx * wy * wz;
}

float normalization_gto(float alpha, uint8_t l, uint8_t m, uint8_t n) {
    double nom =   pow(2.0, 2.0 * (l + m + n) + 3.0 / 2.0) *
                   pow(alpha, (l + m + n) + 3.0 / 2.0);

    double denom = (l < 1 ? 1.0 : double_factorial[(l << 1) - 1] ) *
                   (m < 1 ? 1.0 : double_factorial[(m << 1) - 1] ) *
                   (n < 1 ? 1.0 : double_factorial[(n << 1) - 1] ) *
                   pow(PI, 3.0 / 2.0);

    return sqrt(nom / denom);
}

float overlap_1d(uint8_t l1, uint8_t l2, float x1, float x2, float gamma) {
    float sum = 0.0;

    for(uint8_t i=0; i < (1 + ((l1 + l2) >> 1)); i++) {
        float s = binomial_prefactor(2*i, l1, l2, x1, x2) *
                     (i == 0 ? 1 : double_factorial[(i << 1) - 1] / pow(2 * gamma, i));
        sum += s;
        sprintf(termbuffer, "%i %6.4f", i, s);
        terminal_printtermbuffer();
    }

    return sum;
}

float kinetic_cgf(float *coeff1, float *coeff2, float *alpha1, float *alpha2, float *pos1, float *pos2, uint8_t sz1, uint8_t sz2) {
    float s = 0.0;
    for(uint8_t i=0; i<sz1; i++) {
        for(uint8_t j=0; j<sz2; j++) {
            s += coeff1[i] * coeff2[j] * kinetic(alpha1[i], 0, 0, 0, pos1[0], pos1[1], pos1[2], 
                                                 alpha2[j], 0, 0, 0, pos2[0], pos2[1], pos2[2]);
        }
    }

    return s;
}

float kinetic(float alpha1, uint8_t l1, uint8_t m1, uint8_t n1, float ax, float ay, float az,
              float alpha2, uint8_t l2, uint8_t m2, uint8_t n2, float bx, float by, float bz) {
    float term0 = alpha2 *
                   (2.0 * ( l2 + m2 + n2 ) + 3.0 ) *
                   overlap(alpha1, l1, m1, n1, ax, ay, az, alpha2, l2, m2, n2, bx, by, bz);

    float term1 = -2.0 * pow(alpha2, 2.0) * (
        overlap(alpha1, l1, m1, n1, ax, ay, az, alpha2, l2+2, m2, n2, bx, bx, bz) +
        overlap(alpha1, l1, m1, n1, ax, ay, az, alpha2, l2, m2+2, n2, bx, by, bz) +
        overlap(alpha1, l1, m1, n1, ax, ay, az, alpha2, l2, m2, n2+2, bx, by, bz)
    );

    float term2 = 0.0;
    if(l2 > 1) {
        term2 += -0.5 * l2 * (l2 - 1) * overlap(alpha1, l1, m1, n1, ax, ay, az, alpha2, l2-2, m2, n2, bx, by, bz);
    }
    if(m2 > 1) {
        term2 += -0.5 * m2 * (m2 - 1) * overlap(alpha1, l1, m1, n1, ax, ay, az, alpha2, l2, m2-2, n2, bx, by, bz);
    }
    if(n2 > 1) {
        term2 += -0.5 * n2 * (n2 - 1) * overlap(alpha1, l1, m1, n1, ax, ay, az, alpha2, l2, m2, n2-2, bx, by, bz);
    }

    return term0 + term1 + term2;
}

float binomial_prefactor(int8_t s, int8_t ia, int8_t ib, float xpa, float xpb) {
    float sum = 0.0;

    // sprintf(termbuffer, "Input: %i %i %i %f %f", s, ia, ib, xpa, xpb);
    // terminal_printtermbuffer();

    // sprintf(termbuffer, "s+1: %i", s+1);
    // terminal_printtermbuffer();

    for (int8_t t=0; t<(int8_t)(s+1); t++) {
        // sprintf(termbuffer, "s-ia: %i, s-ia<=t, %i, t, %i", s-ia, (s-ia <= t) ? 1 : 0, t);
        // terminal_printtermbuffer();
        if((s-ia+10 <= t+10) && (t+10 <= ib+10)) { // ultra ugly hack to make this work
            sum += binomial(ia, s-t) *
                   binomial(ib, t) *
                   (ia - s + t) == 0 ? 1.0 : pow(xpa, ia - s + t) *
                   (ib - t) == 0 ? 1.0 : pow(xpb, ib - t);

            // sprintf(termbuffer, "B1: %f", binomial(ia, s-t));
            // terminal_printtermbuffer();
            // sprintf(termbuffer, "B2: %f", binomial(ib, t));
            // terminal_printtermbuffer();
        }
    }

    return sum;
}

float binomial(int8_t a, int8_t b) {
    if( (a < 0) || (b < 0) || (a-b < 0) ) {
        return 1.0;
    }

    return factorial[a] / factorial[b] * factorial[a-b];
}
