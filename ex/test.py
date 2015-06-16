# -*- coding: utf-8 -*-
# enable true divison
#from __future__ import division
import sympy as sp
import rst_symbtools.symb_tools as st

# gehört zu 5.1.3


x1, x2, x3 = sp.symbols("x1, x2, x3")
vec_x = sp.Matrix([x1, x2, x3])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3 = vec_xdot

#F_eq = sp.Matrix([ xdot1 - x3*x4, xdot2 - x4, xdot3 - x5 ])

P1i = sp.Matrix([ [0,1/x3,0],[0,0,1] ])
P0i = sp.Matrix([ [1/x3,-(2*xdot3)/(x3**2),0],[x1/x3, -(2*x1*xdot3)/(x3**2)-(x2/x3)+(2*xdot1)/(x3), 1] ])
