# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symbtools as st

vec_x = st.symb_vector('x1:6')
vec_xdot = st.time_deriv(vec_x, vec_x)

st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)

F_eq = sp.Matrix([ xdot1 - x3*x4, xdot2 - x4, xdot3 - x5 ])
