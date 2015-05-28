# -*- coding: utf-8 -*-
import sympy as sp
import symb_tools as st

a3 = sp.symbols("a3")

x1, x2, x3 = sp.symbols("x1, x2, x3")
vec_x = sp.Matrix([x1, x2, x3])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3 = vec_xdot

F_eq = sp.Matrix([ xdot3 - a3*x1*x2 ])

#P1i = sp.Matrix([ [0,0,1] ])
#P0i = sp.Matrix([ [-a3*x2,-a3*x1,0] ])
