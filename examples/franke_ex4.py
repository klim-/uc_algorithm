# -*- coding: utf-8 -*-
# enable true divison

import sympy as sp
import symbtools as st

x1, x2, x3 = sp.symbols("x1, x2, x3")
vec_x = sp.Matrix([x1, x2, x3])

vec_xdot = st.time_deriv(vec_x, vec_x)
xdot1, xdot2, xdot3 = vec_xdot

F_eq = sp.Matrix([ xdot3 - sp.atan(xdot2/xdot1) ])

#P1i = sp.Matrix([ [(xdot2)/(xdot1**2+xdot2**2), -(xdot1)/(xdot1**2+xdot2**2), 1] ])
#P0i = sp.Matrix([ [0,0,0] ])
