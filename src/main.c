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

#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <z80.h>

#include "config.h"
#include "constants.h"
#include "memory.h"
#include "util.h"
#include "integrals.h"
#include "stack.h"
#include "terminal.h"

// // set printf io
//#pragma printf = classic
#pragma printf "%s %c %i %f"

// forward declarations
void init(void);

int main(void) {
    init();

    float coeff[2][3] = {{0.154329,0.535328,0.444635}, {0.154329,0.535328,0.444635}};
    float alpha[2][3] = {{3.425251,0.623914,0.168855}, {3.425251,0.623914,0.168855}};
    float coeffn[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};    // normalized coefficients
    float pos[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 1.4}};

    sprintf(termbuffer, "coeff: %6.4f %6.4f %6.4f", coeff[0][0], coeff[0][1], coeff[0][2]);
    terminal_printtermbuffer();
    sprintf(termbuffer, "alpha: %6.4f %6.4f %6.4f", alpha[0][0], alpha[0][1], alpha[0][2]);
    terminal_printtermbuffer();

    // normalize coefficients
    for(uint8_t i=0; i<2; i++) {
        for(uint8_t j=0; j<3; j++) {
            coeffn[i][j] = coeff[i][j] * normalization_gto(alpha[i][j], 0, 0, 0);
        }
    }

    float val = 0.0;

    val = overlap(3.425251, 2, 0, 0, 0, 0, 0,
                  3.425251, 2, 0, 0, 0, 0, 0);
    sprintf(termbuffer, "%6.4f", val);
    terminal_printtermbuffer();

    val = overlap(3.425251, 1, 0, 0, 0, 0, 0,
                  3.425251, 1, 0, 0, 0, 0, 0);
    sprintf(termbuffer, "%6.4f", val);
    terminal_printtermbuffer();

    val = overlap(3.425251, 0, 0, 0, 0, 0, 0,
                  3.425251, 0, 0, 0, 0, 0, 0);
    sprintf(termbuffer, "%6.4f", val);
    terminal_printtermbuffer();

    // val = overlap_1d(0, 0, 1, 0, 1);
    // sprintf(termbuffer, "ov1d(0,0,1,0,1) = %6.4f", val);
    // terminal_printtermbuffer();
    // val = overlap_1d(1, 0, 1, 0, 1);
    // sprintf(termbuffer, "ov1d(1,0,1,0,1) = %6.4f", val);
    // terminal_printtermbuffer();
    // val = overlap_1d(2, 0, 1, 0, 1);
    // sprintf(termbuffer, "ov1d(2,0,1,0,1) = %6.4f", val);
    // terminal_printtermbuffer();
    // val = binomial_prefactor(0,2,0,1,0);
    // sprintf(termbuffer, "%6.4f", val);
    // terminal_printtermbuffer();

    //binomial_prefactor(int8_t s, int8_t ia, int8_t ib, float xpa, float xpb)

    print("Calculating overlap integrals:");

    // // calculate overlap integrals
    // for(uint8_t i=0; i<2; i++) {
    //     for(uint8_t j=0; j<2; j++) {
    //         float s = overlap_cgf(coeffn[i], coeffn[j], alpha[i], alpha[j], pos[i], pos[j], 3, 3);
    //         sprintf(termbuffer, "<%i|%i> = %6.4f", i+1, j+1, s);
    //         terminal_printtermbuffer();
    //     }
    // }

    // print("Calculating kinetic integrals:");

    // // calculate kinetic integrals
    // for(uint8_t i=0; i<2; i++) {
    //     for(uint8_t j=0; j<2; j++) {
    //         float s = kinetic_cgf(coeffn[i], coeffn[j], alpha[i], alpha[j], pos[i], pos[j], 3, 3);
    //         sprintf(termbuffer, "<%i|T|%i> = %6.4f", i+1, j+1, s);
    //         terminal_printtermbuffer();
    //     }
    // }

    // sprintf(termbuffer, "%cAll done!", COL_GREEN);
    // terminal_printtermbuffer();

    return 0;
}

void init(void) {
    clear_screen();
    terminal_init(3, 20);
    sprintf(&vidmem[0], "%c%cP2000T Hartree-Fock", TEXT_DOUBLE, COL_CYAN);
    sprintf(&vidmem[0x50*22], "Version: %s", __VERSION__);
    sprintf(&vidmem[0x50*23], "Compiled at: %s / %s", __DATE__, __TIME__);
}
