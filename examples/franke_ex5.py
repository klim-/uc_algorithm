# -*- coding: utf-8 -*-
from __future__ import division
import sympy as sp
import symbtools as st

x1, x2, x3 = sp.symbols("x1, x2, x3")
vec_x = sp.Matrix([x1, x2, x3])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3 = vec_xdot

F_eq = sp.Matrix([ xdot3 + x1*xdot2 - x2*xdot1 ])

#P1i = sp.Matrix([ [x2,-x1,1] ])
#P0i = sp.Matrix([ [-xdot2,xdot1,0] ])
