# -*- coding: utf-8 -*-
import sympy as sp
import symbtools as st
from sympy import sin, cos, tan

# Differentialgleichung 2. Ordnung:
#~ sp.Matrix([
        #~ [m2*q2*(g*sin(p1) + pddot1*q2 + 2*pdot1*qdot2 + qddot1*cos(p1))]
    #~ ])


vec_x = st.symb_vector('x1:7')
vec_xdot = st.time_deriv(vec_x, vec_x)
st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)

g = sp.symbols("g")

diff_symbols = sp.Matrix([])

F_eq = sp.Matrix([
        [ xdot1 - x4 ],
        [ xdot2 - x5 ],
        [ xdot3 - x6 ],
        [ g*sin(x1) + xdot4*x3 + 2*x4*x6 + xdot5*cos(x1)  ]])
