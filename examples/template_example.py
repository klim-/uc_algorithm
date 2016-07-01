# -*- coding: utf-8 -*-
import sympy as sp
import symbtools as st
from sympy import sin, cos, tan

# Number of state variables
n = 6 

vec_x = st.symb_vector('x1:%i' % (n+1))
vec_xdot = st.time_deriv(vec_x, vec_x)
st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)

# Additional symbols
g, l = sp.symbols("g, l")

# Time-dependent symbols
diff_symbols = sp.Matrix([l])

# Nonlinear system in state space representation 0 = F_eq(vec_x, vec_xdot)
F_eq = sp.Matrix([
        [ xdot1 - x4 ],
        [ xdot2 - x5 ],
        [ xdot3 - x6 ],
        [ g*sin(x1) + xdot4*x3 + 2*x4*x6 + xdot5*cos(x1)  ]])
