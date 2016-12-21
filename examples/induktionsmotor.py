# -*- coding: utf-8 -*-
# enable true divison

import sympy as sp
from sympy import sin,cos,tan
import symbtools as st

np, mu, tau, J, M, eta= sp.symbols("np, mu, tau, J, M, eta")

vec_x = st.symb_vector('x1:6')
vec_xdot = st.time_deriv(vec_x, vec_x)

st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)


F_eq = sp.Matrix([
        [xdot1 - mu*x1*x2 + (tau/J)],
        [xdot2 - eta*(M*x4 - x2)],
        [xdot3 - np*x1 - eta*M*(x5/x2)]
    ])

# Ergebnis:
#~ Q-matrix = sp.Matrix([
        
    #~ ])
