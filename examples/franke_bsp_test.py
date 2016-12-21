# -*- coding: utf-8 -*-
# enable true divison

import sympy as sp
import symbtools as st


x1 = st.ExtendedSymbol("x1")
x2 = st.ExtendedSymbol("x2")
x3 = st.ExtendedSymbol("x3")
vec_x = sp.Matrix([x1, x2, x3])

vec_xdot = st.time_deriv(vec_x, vec_x)
xdot1, xdot2, xdot3 = vec_xdot

F_eq = sp.Matrix([ xdot2-xdot1*xdot3 ])


#P1i = sp.Matrix([ [1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0] ])
#P0i = sp.Matrix([ [0,0,-x4,-x3,0],[0,0,0,-1,0],[0,0,0,0,-1] ])
