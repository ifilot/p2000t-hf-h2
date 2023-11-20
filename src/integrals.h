#include <math.h>
#include <stdio.h>
#include "memory.h"

#define PI 3.141527
#pragma printf "%i %f"

float factorial(uint8_t i);
float double_factorial(uint8_t i);
float overlap_1d(uint8_t l1, uint8_t l2, float x1, float x2, float gamma);
float overlap(float alpha1, uint8_t l1, uint8_t m1, uint8_t n1, float ax, float ay, float az,
              float alpha2, uint8_t l2, uint8_t m2, uint8_t n2, float bx, float by, float bz);
float binomial(uint8_t a, uint8_t b);
float binomial_prefactor(uint8_t s, uint8_t ia, uint8_t ib, float xpa, float xpb);
