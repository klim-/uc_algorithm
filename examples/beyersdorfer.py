# -*- coding: utf-8 -*-
# enable true divison

import sympy as sp
from sympy import sin,cos,tan
import symbtools as st

L1, L2 = sp.symbols("L1,L2", commutative=False)

vec_x = st.symb_vector('x1:7', commutative=False)
vec_xdot = st.time_deriv(vec_x, vec_x)

st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)


# zur kontrolle
#~ P1 = sp.Matrix([
        #~ [                                -tan(x4), 1, 0,  0, 0, 0],
        #~ [             -tan(x3)*(L1*cos(x4))**(-1), 0, 0,  1, 0, 0],
        #~ [-sin(x5 + x6)*(L2*cos(x5)*cos(x4))**(-1), 0, 0, -1, 0, 1]
    #~ ])

#~ P0 = sp.Matrix([
        #~ [0, 0,                                          0,                                                                    -xdot1 - xdot1*tan(x4)**2,                                                                                                                                            0,                                              0],
        #~ [0, 0, -xdot1*(1 + tan(x3)**2)*(L1*cos(x4))**(-1),                              -xdot1*tan(x3)*(L1*cos(x4))**(-1)*L1*sin(x4)*(L1*cos(x4))**(-1),                                                                                                                                            0,                                              0],
        #~ [0, 0,                                          0, -xdot1*sin(x5 + x6)*(L2*cos(x5)*cos(x4))**(-1)*L2*cos(x5)*sin(x4)*(L2*cos(x5)*cos(x4))**(-1), -xdot1*sin(x5 + x6)*(L2*cos(x5)*cos(x4))**(-1)*L2*sin(x5)*cos(x4)*(L2*cos(x5)*cos(x4))**(-1) - xdot1*cos(x5 + x6)*(L2*cos(x5)*cos(x4))**(-1), -xdot1*cos(x5 + x6)*(L2*cos(x5)*cos(x4))**(-1)]
#~ ])

F_eq = sp.Matrix([
        [xdot2 - xdot1*sp.tan(x4)],
        [xdot4 - (xdot1*sp.tan(x3))/(L1*sp.cos(x4))],
        [xdot6 - xdot4 -(xdot1*sp.sin(x5+x6))/(L2*sp.cos(x5)*sp.cos(x4)) ]
    ])

# Ergebnis:
#~ Q-matrix = sp.Matrix([
        #~ [                                -tan(x4), 1, 0,  0, 0, 0],
        #~ [             -tan(x3)*(L1*cos(x4))**(-1), 0, 0,  1, 0, 0],
        #~ [-sin(x5 + x6)*(L2*cos(x5)*cos(x4))**(-1), 0, 0, -1, 0, 1]
    #~ ])
