# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symb_tools as st

L, C, U, R = sp.symbols("L, C, U, R")

x1, x2 = sp.symbols("x1, x2")
vec_x = sp.Matrix([x1, x2])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2 = vec_xdot


F_eq = sp.Matrix([ L*x1*xdot1 + C*x2*xdot2 - x1*U + (1/R)*x2**2 ])

#P1i = sp.Matrix([ [L*x1, C*x2] ])
#P0i = sp.Matrix([ [L*xdot1-U, C*xdot2+(2/R)*x2] ])
