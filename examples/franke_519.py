# -*- coding: utf-8 -*-
import sympy as sp
import rst_symbtools.symb_tools as st

a0, a1, a2 = sp.symbols("a0, a1, a2")
b0, b1, b2 = sp.symbols("b0, b1, b2")
c0, c1, c2 = sp.symbols("c0, c1, c2")

diff_symbols = sp.Matrix([a0, a1, a2, b0, b1, b2, c0, c1, c2])

x1, x2, x3, x4, x5, x6 = sp.symbols("x1, x2, x3, x4, x5, x6")
vec_x = sp.Matrix([x1, x2, x3, x4, x5, x6])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4, xdot5, xdot6 = vec_xdot

F_eq = sp.Matrix([
        [ xdot1 - x4 ],
        [ xdot2 - x5 ],
        [ xdot3 - x6 ],
        [ a2*xdot4 + b2*xdot5 + c2*xdot6 + b1*x5 + c1*x6 + b0*x2 + c0*x3 ]])
        #~ [ a2*xdot4 + b2*xdot5 + c2*xdot6 + b1*x5 + c1*x6 + b0*x2 + c0*x3 + a0*x1 + a1*x4]])
