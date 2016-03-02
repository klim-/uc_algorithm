# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symbtools as st

T_, K_L, K_V, A_B, A_sp, g, K_S, d, m = sp.symbols("T, K_L, K_V, A_B, A_sp, g, K_S, d, m")

x1, x2, x3, x4 = sp.symbols("x1:5")
vec_x = sp.Matrix([x1, x2, x3, x4])

vec_xdot = st.time_deriv(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4 = vec_xdot

F_eq = sp.Matrix([
[ xdot1 - x2 ],
[ xdot3 - x4 ],
[ xdot4 - (K_L/m)*( (K_V*x1 - A_B*x4)/(A_sp) )**2 + g ]])
