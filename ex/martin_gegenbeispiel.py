# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symb_tools as st

# set up the system
s = sp.symbols("s")

x1, x2, x3, x4, x5, x6, x7 = sp.symbols("x1, x2, x3, x4, x5, x6, x7")
vec_x = sp.Matrix([x1, x2, x3, x4, x5, x6, x7])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4, xdot5, xdot6, xdot7 = vec_xdot


F_eq = sp.Matrix([
[ x3 + x4*xdot1 - xdot2 ],
[ -x6 + x7*x4*xdot1 - xdot5 ],
[ -x7 - x5*xdot1 - x7*xdot3 + (x4 + xdot1)*x4*xdot1 - xdot6 ],
[ x4 + xdot1 - xdot7 ]])
