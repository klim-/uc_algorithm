# -*- coding: utf-8 -*-
# enable true divison

import sympy as sp
import symbtools as st

x1, x2 = sp.symbols("x1, x2")
vec_x = sp.Matrix([x1,x2])

vec_xdot = st.time_deriv(vec_x, vec_x)
xdot1, xdot2 = vec_xdot

F_eq = sp.Matrix([
[xdot1 + x1 + xdot2*x1**2]
])


#P1i = sp.Matrix([ [1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0] ])
#P0i = sp.Matrix([ [0,0,-x4,-x3,0],[0,0,0,-1,0],[0,0,0,0,-1] ])
