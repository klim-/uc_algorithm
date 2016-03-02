# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symbtools as st

Cmms, Ciis, Csis, Csms, Mm, eta, tau, V, alpha1, alpha4 =\
        sp.symbols("Cmms, Ciis, Csis, Csms, Mm, eta, tau, V, alpha1, alpha4")

Rm, ki, theta, f6 = sp.symbols("Rm, ki, theta, f6")
diff_symbols = sp.Matrix([Rm, ki, theta, f6])


vec_x = st.symb_vector('x1:7')
vec_xdot = st.time_deriv(vec_x, vec_x)

st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)



delta = (1 - eta*(x4/(x4 + Mm*x1)))

#~ F_eq = sp.Matrix([
        #~ [xdot1 - (1/tau)*(Cmms - delta*x1)],
        #~ [xdot2 + ki*x2 - (Ciis/Csis)*(xdot3 - (Csms/tau) + delta*(x3/tau)) + delta*(x2/tau)],
        #~ [xdot4 + Mm*Rm + delta*(x4/tau)],
        #~ [xdot5 - theta - alpha1*x6]
    #~ ])

b1,b2,b3,b4 = sp.symbols("b1,b2,b3,b4")
diff_symbols = sp.Matrix([b1,b2,b3,b4])
b5,b6,b7,b8 = sp.symbols("b5,b6,b7,b8")
diff_symbols = sp.Matrix([b5,b6,b7,b8])

F_eq = sp.Matrix([
        [xdot1 - x4],
        [xdot2 - x3],
        [xdot4 - (b1*xdot1 + b2*xdot2 + b3*xdot4 + b4*xdot5 + b5*x1 + b6*x2 + b7*x4 + b8*x5)],
        [xdot5 - x6]
    ])

