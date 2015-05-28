# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symb_tools as st

x1, x2 = st.symbols("x1, x2")
vec_x = sp.Matrix([x1, x2])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2 = vec_xdot

# mit zwei multipliziert damit kein bruch in gleichung
F_eq = sp.Matrix([ 
[ -2*xdot2 + xdot1**2 ]])

#P1i = sp.Matrix([ [-xdot1,1] ])
#P0i = sp.Matrix([ [0,0] ])
