# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symbtools as st

x1, x2, x3 = sp.symbols("x1, x2, x3")
vec_x = sp.Matrix([x1, x2, x3])

vec_xdot = st.time_deriv(vec_x, vec_x)
xdot1, xdot2, xdot3 = vec_xdot

F_eq = sp.Matrix([ xdot1*sp.sin(x3) - xdot2*sp.cos(x3) ])

#P1i = sp.Matrix([ [sp.sin(x3), -sp.cos(x3), 0] ])
#P0i = sp.Matrix([ [0, 0, xdot1*sp.cos(x3)+xdot2*sp.sin(x3)] ])
