# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symb_tools as st

x1, x2, x3, x4 = sp.symbols("x1, x2, x3, x4")
vec_x = sp.Matrix([x1, x2, x3, x4])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4 = vec_xdot

F_eq = sp.Matrix([ xdot1*sp.sin(x3) - xdot2*sp.cos(x3),
                   sp.sqrt(xdot1**2 + xdot2**2) - xdot4 ])

#P1i = sp.Matrix([ [sp.sin(x3), -sp.cos(x3), 0] ])
#P0i = sp.Matrix([ [0, 0, xdot1*sp.cos(x3)+xdot2*sp.sin(x3)] ])
