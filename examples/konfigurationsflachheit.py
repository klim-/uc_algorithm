# -*- coding: utf-8 -*-
import sympy as sp
import symb_tools as st

# set up the system
s = sp.symbols("s")
c1, d1, m2 = sp.symbols("c1, d1, m2")

x1, x2, x3, x4 = sp.symbols("x1, x2, x3, x4")
vec_x = sp.Matrix([x1, x2, x3, x4])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4 = vec_xdot


F_eq = sp.Matrix([ xdot1-x3, xdot2-x4, xdot4+(d1/m2)*x3 - (d1/m2)*x4 + (c1/m2)*x1 - (c1/m2)*x2 ])

#P1i = sp.Matrix([ [1,0,0,0],[0,1,0,0],[0,0,0,1] ])
#P0i = sp.Matrix([ [0,0,-1,0],[0,0,0,-1],[c1/m2,-c1/m2,d1/m2,-d1/m2] ])



