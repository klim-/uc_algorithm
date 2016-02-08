# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symb_tools as st

# set up the system
s = sp.symbols("s")

vec_x = st.symb_vector('x1:8', commutative=True)
vec_xdot = st.perform_time_derivative(vec_x, vec_x)
vec_xddot = st.perform_time_derivative(vec_xdot, vec_xdot)

st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)
st.make_global(vec_xddot, 1)


F_eq = sp.Matrix([
[ x3 + x4*xdot1 - xdot2 ],
[ -x6 + x7*x4*xdot1 - xdot5 ],
[ -x7 - x5*xdot1 - x7*xdot3 + (x4 + xdot1)*x4*xdot1 - xdot6 ],
[ x4 + xdot1 - xdot7 ]])


Q = sp.Matrix([
[x4*x7*xdot4*(x4*x7*xdot1 + 2*x4*xddot1 + 2*x4*xdot4 + 2*xdot1*xdot4 + xdot3 - xdot4*xdot7 - xdot5 + 1) + x4*xdot4*(x4**2*xdot7 - x4*x7**2*xdot1 - 2*x4*x7*xddot1 - 2*x4*x7*xdot4 + x4*xdot1*xdot7 + x4*xdot7**2 - x5*xdot7 - 2*x7*xdot1*xdot4 - x7*xdot3 + x7*xdot4*xdot7 + x7*xdot5 - x7) + xdot1*xdot4*(2*x4 + xdot1 - xdot7)*(x4*xdot7 + x4*(x4 + xdot1) - x5) - xdot4*(x4*xdot1 + x4*(x4 + xdot1) - x5)*(x4*xdot7 + x4*(x4 + xdot1) - x5), -xdot4*(x4**2*xdot7 - x4*x7**2*xdot1 - 2*x4*x7*xddot1 - 2*x4*x7*xdot4 + x4*xdot1*xdot7 + x4*xdot7**2 - x5*xdot7 - 2*x7*xdot1*xdot4 - x7*xdot3 + x7*xdot4*xdot7 + x7*xdot5 - x7), x7*xdot4*(x4*xdot7 + x4*(x4 + xdot1) - x5), 0, -xdot4*(x4*x7*xdot1 + 2*x4*xddot1 + 2*x4*xdot4 + 2*xdot1*xdot4 + xdot3 - xdot4*xdot7 - xdot5 + 1), xdot4*(x4*xdot7 + x4*(x4 + xdot1) - x5), -xdot1*xdot4*(2*x4 + xdot1 - xdot7)*(x4*xdot7 + x4*(x4 + xdot1) - x5)],
[                                                                                                                                                                                                                                                                                                                                                                                                              x4*xdot1 + x4*(x4 + xdot1) - x5,                                                                                                                                                                               0,                                        -x7, 0,                                                                                                 0,                                      -1,                                                                     0],
[                                                                                                                                                                                                                                                                                                                                                                                                                                            1,                                                                                                                                                                               0,                                          0, 0,                                                                                                 0,                                       0,                                                                    -1]])

