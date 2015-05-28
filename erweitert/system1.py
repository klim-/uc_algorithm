# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import rst_symbtools.symb_tools as st

m1 ,m2, d1, d0, c1 = sp.symbols("m1 ,m2, d1, d0, c1")

x1 = st.ExtendedSymbol("x1")
x2 = st.ExtendedSymbol("x2")
vec_x = sp.Matrix([x1, x2])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2 = vec_xdot


P2i = sp.Matrix([
[0,m2]])

P1i = sp.Matrix([
[d1,-d1]])

P0i = sp.Matrix([
[c1,-c1]])
