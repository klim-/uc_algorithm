# -*- coding: utf-8 -*-
# enable true divison

import sympy as sp
import symbtools as st

vec_x = st.symb_vector('x1:6')
vec_xdot = st.time_deriv(vec_x, vec_x)

st.make_global(vec_x)
st.make_global(vec_xdot)

F_eq = sp.Matrix([
        [-x3*x4 + xdot1],
        [-x4 + xdot2],
        [-x5 + xdot3]])
