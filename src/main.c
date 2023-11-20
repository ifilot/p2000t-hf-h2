#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "constants.h"
#include "memory.h"
#include "util.h"
#include "integrals.h"
#include "stack.h"
#include "z80.h"

#define LOWMEM 0x6200
#define BANKMEM 0xE000
#define BANK_BYTES 0x2000
#define STACK 0x9F00

// forward declarations
void init(void);

int main(void) {
    init();

    float coeff[3] = {0.154329,0.535328,0.444635};
    float alpha[3] = {3.425251,0.623914,0.168855};
    float pos1[3] = {0.0, 0.0, 0.0};
    float pos2[3] = {0.0, 0.0, 0.0};

    sprintf(&vidmem[0x00], "coeff: %6.4f %6.4f %6.4f", coeff[0], coeff[1], coeff[2]);
    sprintf(&vidmem[0x50], "alpha: %6.4f %6.4f %6.4f", alpha[0], alpha[1], alpha[2]);
    sprintf(&vidmem[0x50*2], "bnpf: %6.4f", binomial_prefactor(0,0,0,0.0,0.0));
    sprintf(&vidmem[0x50*3], "bn: %6.4f", binomial(0,0));
    sprintf(&vidmem[0x50*4], "pow: %6.4f", pow(0.0, 0));

    float o00 = overlap(alpha[0], 0, 0, 0, pos1[0], pos1[1], pos1[2],
                        alpha[0], 0, 0, 0, pos1[0], pos1[1], pos1[2]);
    sprintf(&vidmem[0x50*5], "overlap gto1: %6.4f", o00);

    float res = 0.0;
    for(uint8_t i=0; i<3; i++) {
        for(uint8_t j=0; j<3; j++) {
            float s = overlap(alpha[i], 0, 0, 0, pos1[0], pos1[1], pos1[2],
                              alpha[j], 0, 0, 0, pos2[0], pos2[1], pos2[2]);
            res += coeff[i] * coeff[j] * s;
            sprintf(&vidmem[0x50*12], "%i %i %6.4f", i, j, s);
        }
    }
    sprintf(&vidmem[0x50*13], "Result of overlap integral: %6.4f", res);

    return 0;
}

void init(void) {
    clearline(0);
    sprintf(&vidmem[0x50*22], "Version: %s", __VERSION__);
    sprintf(&vidmem[0x50*23], "Compiled at: %s / %s", __DATE__, __TIME__);
}
