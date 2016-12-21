# -*- coding: utf-8 -*-
# enable true divison

import sympy as sp
import symbtools as st

C,L,Ld,Rd,Ll,Rl,Ul,U = sp.symbols("C,L,Ld,Rd,Ll,Rl,Ul,U", commutative=True)

vec_x = st.symb_vector('x1:13', commutative=True)
vec_xdot = st.time_deriv(vec_x, vec_x)

st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)


F_eq = sp.Matrix([
        #~ [L*xdot5 + x9*x1 + x10*x2 + L*xdot6 - U + Ld*xdot5 + Ld*xdot7 - Rd*x5 - Rd*x7],
        [ xdot8 - (1-Rd/Rl)*(x5+x7) ],
        [L*xdot7 + x11*x3 + x12*x4 + L*xdot8 - U + Ld*xdot5 + Ld*xdot7 + Rd*x5 + Rd*x7],
        [L*xdot5 - L*xdot7 - x11*x3 + x9*x1 + Ll*xdot5 - Ll*xdot6 + Rl*x5 + Rl*x6 + Ul],
        [L*xdot6 - L*xdot8 + x10*x2 - x12*x4 - Ll*xdot5 + Ll*xdot6 - Rl*x5 + Rl*x6 - Ul],
        [C*xdot1 - x5*x9],
        [C*xdot2 - x6*x10],
        [C*xdot3 - x7*x11],
        [C*xdot4 - x8*x12],
        [ xdot6 - (Rd/Rl)*(xdot5 + xdot7) ]
        #~ [x5 + x7 - x6 - x8]
    ])

# algebraische Substitutionsgleichungen mit Ableitungen:
x6subs = -Rd*(x5 + x7)/Rl
xdot6subs = -Rd*(xdot5 + xdot7)/Rl
x8subs =  (x5 + x7)*(Rd/Rl + 1)
xdot8subs = (xdot5 + xdot7)*(Rd/Rl + 1)

F_eq_new = F_eq.subs(x8,x8subs).subs(xdot8,xdot8subs).subs(x6,x6subs).subs(xdot6,xdot6subs)

P10 = F_eq_new.jacobian(vec_xdot)

