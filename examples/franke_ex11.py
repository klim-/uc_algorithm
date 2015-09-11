# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symb_tools as st

a = sp.symbols("a")

x1, x2, x3, x4, x5 = sp.symbols("x1, x2, x3, x4, x5")
vec_x = sp.Matrix([x1, x2, x3, x4, x5])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4, xdot5 = vec_xdot

F_eq = sp.Matrix([
[ xdot1*sp.cos(x5) + xdot2*sp.sin(x5) + a*(xdot5*sp.cos(x3) + xdot4) ],
[ -xdot1*sp.sin(x5) + xdot2*sp.cos(x5) + a*xdot3*sp.sin(x3) ] ])

#P1i = sp.Matrix([ [sp.cos(x5),sp.sin(x5),0,a,a*sp.cos(x3)],[-sp.sin(x5),sp.cos(x5),a*sp.sin(x3),0,0] ])
#P0i = sp.Matrix([ [0,0,-a*xdot5*sp.sin(x3),0,xdot2*sp.cos(x5)-xdot1*sp.sin(x5)],[0,0,a*xdot3*sp.cos(x3),0,xdot1*sp.cos(x5)+xdot2*sp.sin(x5)] ])
