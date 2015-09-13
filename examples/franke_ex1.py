# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symb_tools as st

x1, x2, x3, x4, x5 = sp.symbols("x1, x2, x3, x4, x5")
vec_x = sp.Matrix([x1, x2, x3, x4, x5])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4, xdot5 = vec_xdot

F_eq = sp.Matrix([ xdot1 - x3*x4, xdot2 - x4, xdot3 - x5 ])
