# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symbtools as st

k, m, R_m, R_e, b, L_x, L = sp.symbols("k, m, R_m, R_e, b, L_x, L")

x1, x2, x3 = xx = sp.symbols("x1:4")
vec_x = sp.Matrix([x1, x2, x3])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3 = vec_xdot

F_eq = sp.Matrix([
[ xdot1 - x2 ],
[ xdot2 + (k/m)*x1 + (R_m/m)*x2 - (b/m)*x3 + (L_x/(2*m))*x3**2 ]])
