# -*- coding: utf-8 -*-
# enable true divison
#from __future__ import division
import sympy as sp
import rst_symbtools.symb_tools as st



x1, x2 = sp.symbols("x1, x2")
vec_x = sp.Matrix([x1, x2])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2 = vec_xdot

#F_eq = sp.Matrix([-5*x1 + x1**3 - (2/3)*xdot1 - xdot2])
F_eq = sp.Matrix([-5*x1**2 + x1 - 2*xdot1-xdot2])

