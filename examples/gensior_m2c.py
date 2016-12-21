# -*- coding: utf-8 -*-
# enable true divison

import sympy as sp
from sympy import sin,cos,tan
import symbtools as st

C,L,Ld,Rd,Ll,Rl,Ul,U = sp.symbols("C,L,Ld,Rd,Ll,Rl,Ul,U", commutative=True)

vec_x = st.symb_vector('x1:13', commutative=True)
vec_xdot = st.time_deriv(vec_x, vec_x)

st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)


# zur kontrolle
#~ P1 = sp.Matrix([
        #~ Matrix([
        #~ [0, 0, 0, 0, 2*L + 2*Ld, 0, 2*L + 2*Ld, 0, 0, 0, 0],
        #~ [0, 0, 0, 0,          0, 0,         -L, 0, 0, 0, 0],
        #~ [C, 0, 0, 0,          0, 0,          0, 0, 0, 0, 0],
        #~ [0, C, 0, 0,          0, 0,          0, 0, 0, 0, 0],
        #~ [0, 0, C, 0,          0, 0,          0, 0, 0, 0, 0],
        #~ [0, 0, 0, C,          0, 0,          0, 0, 0, 0, 0]])

    #~ ])

#~ P0 = sp.Matrix([
        #~ [x8,  x9,  x10,  x11, 2*Rd,   0, 2*Rd,  x1,  x2,  x3,            x4],
        #~ [x8, -x9, -x10, -x11,    0,   0,    0,  x1, -x2, -x3,           -x4],
        #~ [ 0,   0,    0,    0,  -x8,   0,    0, -x5,   0,   0,             0],
        #~ [ 0,   0,    0,    0,    0, -x9,    0,   0, -x6,   0,             0],
        #~ [ 0,   0,    0,    0,    0,   0, -x10,   0,   0, -x7,             0],
        #~ [ 0,   0,    0,    0, -x11, x11, -x11,   0,   0,   0, -x5 + x6 - x7]])
    #~ ])

# alte gleichungen komplett:
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

x6subs = -Rd*(x5 + x7)/Rl
xdot6subs = -Rd*(xdot5 + xdot7)/Rl
x8subs =  (x5 + x7)*(Rd/Rl + 1)
xdot8subs = (xdot5 + xdot7)*(Rd/Rl + 1)

F_eq_new = F_eq.subs(x8,x8subs).subs(xdot8,xdot8subs).subs(x6,x6subs).subs(xdot6,xdot6subs)

P10 = F_eq_new.jacobian(vec_xdot)

# x6 und x8 eliminiert:
#~ alt:     (x1,x2,x3,x4,x5,x7,x9,x10,x11,x12)
#~ neu:-->  (x1,x2,x3,x4,x5,x6,x7, x8, x9,x10)
#~ F_eq = sp.Matrix([
        #~ [ x3*x9 + x4*x10 + L*xdot6 + (L*(1+Rd/Rl)+Rd)*(xdot5+xdot6) + Rd*(x5+x6) ],
        #~ [ x3*x9 + x1*x7 + (L+Ll+Rd/Rl)*xdot5 + (Ll*Rd/Rl-L)*xdot6 + (Rl+Rd)*x5 + Rd*x6 + Ul ],
        #~ [ (Ll+Rd/Rl - L - 2*L*Rd/Rl)*xdot5 + (Ll*Rd/Rl - L - 2*Rd/Rl)*xdot6 + x2*x8 - x4*x10 - (Rl+Rd)*x5 - Rd*x6 - Ul ],
        #~ [ C*xdot1 - x5 * x7 ],
        #~ [ C*xdot2 + (Rd/Rl)*(x5+x6)*x8 ],
        #~ [ C*xdot3 - x6*x9 ],
        #~ [C*xdot4 - (1-Rd/Rl)*(x5+x6)*x10 ]
    #~ ])



# Ergebnis:
#~ Q-matrix = sp.Matrix([
        #~ Q-matrix = 
        #~ Matrix([
        #~ [2*C*x1, 0, 0, 0, x5*(2*L + 2*Ld), 0, -L*x5 + x5*(2*L + 2*Ld), 0, 0, 0, 0],
        #~ [     C, 0, 0, 0,               0, 0,                       0, 0, 0, 0, 0],
        #~ [     0, C, 0, 0,               0, 0,                       0, 0, 0, 0, 0],
        #~ [     0, 0, C, 0,               0, 0,                       0, 0, 0, 0, 0],
        #~ [     0, 0, 0, C,               0, 0,                       0, 0, 0, 0, 0]])
    #~ ])

#~ A = sp.Matrix([
        #~ [L+Ld, L, Ld, 0],
        #~ [Ld,0,L+Ld,L],
        #~ [L+Ll,-Ll,-L,0],
        #~ [-Ll,L+Ll,0,-L]])
