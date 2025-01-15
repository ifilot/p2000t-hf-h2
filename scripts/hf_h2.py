# -*- coding: utf-8 -*-

from pyqint import PyQInt, cgf, gto, Molecule
from copy import deepcopy
import numpy as np

# construct integrator object
integrator = PyQInt()

cgf1 = cgf([0.0, 0.0, 0.0])
cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

cgf2 = deepcopy(cgf1)
cgf2.p[2] = 1.4

g1 = gto(1.0, [0,0,0], 3.425251, 0, 0, 0)
print(integrator.overlap_gto(g1, g1))

g2 = gto(1.0, [0,0,0], 3.425251, 1, 0, 0)
print(integrator.overlap_gto(g2, g2))


g2 = gto(1.0, [0,0,0], 3.425251, 2, 0, 0)
print(integrator.overlap_gto(g2, g2))


# S = np.zeros((2,2))
# S[0,0] = integrator.overlap(cgf1, cgf1)
# S[0,1] = S[1,0] = integrator.overlap(cgf1, cgf2)
# S[1,1] = integrator.overlap(cgf2, cgf2)

#print(S)