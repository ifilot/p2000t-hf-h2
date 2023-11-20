#include "integrals.h"

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

    sprintf(&vidmem[0x50*2], "%6.4f %6.4f %6.4f %6.4f", pre, wx, wy, wz);

    return pre * wx * wy * wz;
}

float overlap_1d(uint8_t l1, uint8_t l2, float x1, float x2, float gamma) {
    float sum = 0.0;

    for(uint8_t i=0; i < (1 + floor(0.5 * (l1 + l2))); i++) {
        sum += binomial_prefactor(2*i, l1, l2, x1, x2) *
                     (i == 0 ? 1 : double_factorial(2 * i - 1) ) /
                     pow(2 * gamma, i);
    }

    return sum;
}

float binomial_prefactor(uint8_t s, uint8_t ia, uint8_t ib, float xpa, float xpb) {
    float sum = 0.0;

    for (uint8_t t=0; t < s+1; t++) {
        if ((s-ia <= t) && (t <= ib)) {
            sum += binomial(ia, s-t) *
                   binomial(ib, t) *
                   (ia - s + t) == 0 ? 1.0 : pow(xpa, ia - s + t) *
                   (ib - t) == 0 ? 1.0 : pow(xpb, ib - t);
        }
    }

    return sum;
}

float binomial(uint8_t a, uint8_t b) {
    return factorial(a) / factorial(b) * factorial(a-b);
}

float factorial(uint8_t i) {
    static const float res[] = {1,1,2,6,24,120,720,5040,40320,362880,3628800};
    return res[i];
}

float double_factorial(uint8_t i) {
    static const float res[] = {1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840};
    return res[i];
}
