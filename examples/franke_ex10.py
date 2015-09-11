# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symb_tools as st

e = sp.symbols("epsilon")

x1, x2, x3, x4, x5, x6 = sp.symbols("x1, x2, x3, x4, x5, x6")
vec_x = sp.Matrix([x1, x2, x3, x4, x5, x6])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4, xdot5, xdot6 = vec_xdot


F_eq = sp.Matrix([ [xdot4*sp.cos(x3) - (xdot5+1)*sp.sin(x3) + e*xdot6], [xdot1-x4], [xdot2-x5], [xdot3-x6] ])

#P1i = sp.Matrix([ [0,0,0,sp.cos(x3),-sp.sin(x3),eta],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0] ])
#P0i = sp.Matrix([ [0,0,-xdot4*sp.sin(x3)-(xdot5+1)*sp.cos(x3),0,0,0],[0,0,0,-1,0,0],[0,0,0,0,-1,0],[0,0,0,0,0,-1] ])
