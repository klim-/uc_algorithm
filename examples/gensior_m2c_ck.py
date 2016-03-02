# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
import symbtools as st

C,L,Ld,Rd,Ll,Rl,Ul,U = sp.symbols("C,L,Ld,Rd,Ll,Rl,Ul,U", commutative=True)

xx = vec_x = st.symb_vector('x1:13', commutative=True)
xxd = vec_xdot = st.perform_time_derivative(vec_x, vec_x)

st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)


#~ F_eq = sp.Matrix([
        #~ [L*xdot5 + x9*x1 + x10*x2 + L*xdot6 - U + Ld*xdot5 + Ld*xdot7 - Rd*x5 - Rd*x7],
        #~ #[ xdot8 - (1-Rd/Rl)*(x5+x7) ],
        #~ [L*xdot7 + x11*x3 + x12*x4 + L*xdot8 - U + Ld*xdot5 + Ld*xdot7 + Rd*x5 + Rd*x7],
        #~ [L*xdot5 - L*xdot7 - x11*x3 + x9*x1 + Ll*xdot5 - Ll*xdot6 + Rl*x5 + Rl*x6 + Ul],
        #~ [L*xdot6 - L*xdot8 + x10*x2 - x12*x4 - Ll*xdot5 + Ll*xdot6 - Rl*x5 + Rl*x6 - Ul],
        #~ [C*xdot1 - x5*x9],
        #~ [C*xdot2 - x6*x10],
        #~ [C*xdot3 - x7*x11],
        #~ [C*xdot4 - x8*x12],
        #~ #[ xdot6 - (Rd/Rl)*(xdot5 + xdot7) ],
        #~ [x5 + x7 - x6 - x8]
    #~ ])

#~ P1a = F_eq.jacobian(xxd)
#~ NSa = st.nullspaceMatrix(P1a.T).T

#~ # Algebraische Gleichungen
#~ eq_alg_a = NSa*F_eq

#~ eq_alg_a_dot = st.perform_time_derivative(eq_alg_a, xx)


#~ rest_eq = st.nullspaceMatrix(NSa).T*F_eq

#~ Fb = st.row_stack(rest_eq, eq_alg_a_dot)

#~ st.rnd_number_rank(Fb.jacobian(xxd)) # -> 8
#~ P1b = Fb.jacobian(xxd)


#~ NSb = st.nullspaceMatrix(P1b.T).T

#~ NSb.simplify()


#~ # wieder die algebraischen Gleichungen bestimmen
#~ eq_alg_b = NSb * Fb

#~ # Probe, dass keine Ableitungen mehr auftreten (simplify alleine kriegt es nicht hin)
#~ assert sp.simplify(eq_alg_b.jacobian(xxd)) == 0*eq_alg_b.jacobian(xxd)

#~ # Ableitungen direkt auf 0 setzen (manuelles vereinfachen)
#~ eq_alg_b = eq_alg_b.subz0(xxd)
#~ eq_alg_b_dot = st.perform_time_derivative(eq_alg_b, xx)


#~ # neue Systemgleichungen
#~ Fc = st.row_stack(Fb[1:,:], eq_alg_b_dot)

# ergibt:
#~ Fc = sp.Matrix([
        #~ [                                                                                                                                                                                                                                                                                                                                       2*L*xdot5 + L*xdot6 - L*xdot7 + Ld*xdot5 + Ld*xdot7 + Ll*xdot5 - Ll*xdot6 - Rd*x5 - Rd*x7 + Rl*x5 + Rl*x6 - U + Ul + 2*x1*x9 + x10*x2 - x11*x3],
        #~ [                                                                                                                                                                                                                                                                                                                                       L*xdot5 + 2*L*xdot6 - L*xdot8 + Ld*xdot5 + Ld*xdot7 - Ll*xdot5 + Ll*xdot6 - Rd*x5 - Rd*x7 - Rl*x5 + Rl*x6 - U - Ul + x1*x9 + 2*x10*x2 - x12*x4],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                                      C*xdot1 - x5*x9],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                                     C*xdot2 - x10*x6],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                                     C*xdot3 - x11*x7],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                                     C*xdot4 - x12*x8],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                 2*Rd*xdot5 + 2*Rd*xdot7 + 2*Rl*xdot6],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                        xdot5 - xdot6 + xdot7 - xdot8],
        #~ [Rl*(-L**2 - L*Ld + L*Ll + Ld*Ll)*(-Rd*xdot5 - Rd*xdot7 + Rl*xdot5 + Rl*xdot6 + 2*x1*xdot9 + x10*xdot2 - x11*xdot3 + x2*xdot10 - x3*xdot11 + 2*x9*xdot1) + Rl*(3*L**2 + 3*L*Ld + L*Ll + Ld*Ll)*(-Rd*xdot5 - Rd*xdot7 - Rl*xdot5 + Rl*xdot6 + x1*xdot9 + 2*x10*xdot2 - x12*xdot4 + 2*x2*xdot10 - x4*xdot12 + x9*xdot1) + (4*L**2*Rd + L**2*Rl - L*Ld*Rl + 4*L*Ll*Rd + L*Ll*Rl - Ld*Ll*Rl)*(x1*xdot9 + x10*xdot2 + x11*xdot3 + x12*xdot4 + x2*xdot10 + x3*xdot11 + x4*xdot12 + x9*xdot1)]
    #~ ])

# Konstanten zusammenfassen
# K1 = Rl*(-L**2 - L*Ld + L*Ll + Ld*Ll)
# K2 = Rl*(3*L**2 + 3*L*Ld + L*Ll + Ld*Ll)
# K3 = (4*L**2*Rd + L**2*Rl - L*Ld*Rl + 4*L*Ll*Rd + L*Ll*Rl - Ld*Ll*Rl)
#~ K1, K2, K3 = sp.symbols("K1, K2, K3", commutative=True)
#~ Fc = sp.Matrix([
        #~ [                                                                                                                                                                                                                                                                                                                                       2*L*xdot5 + L*xdot6 - L*xdot7 + Ld*xdot5 + Ld*xdot7 + Ll*xdot5 - Ll*xdot6 - Rd*x5 - Rd*x7 + Rl*x5 + Rl*x6 - U + Ul + 2*x1*x9 + x10*x2 - x11*x3],
        #~ [                                                                                                                                                                                                                                                                                                                                       L*xdot5 + 2*L*xdot6 - L*xdot8 + Ld*xdot5 + Ld*xdot7 - Ll*xdot5 + Ll*xdot6 - Rd*x5 - Rd*x7 - Rl*x5 + Rl*x6 - U - Ul + x1*x9 + 2*x10*x2 - x12*x4],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                                      C*xdot1 - x5*x9],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                                     C*xdot2 - x10*x6],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                                     C*xdot3 - x11*x7],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                                     C*xdot4 - x12*x8],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                 2*Rd*xdot5 + 2*Rd*xdot7 + 2*Rl*xdot6],
        #~ [                                                                                                                                                                                                                                                                                                                                                                                                                                                        xdot5 - xdot6 + xdot7 - xdot8],
        #~ [K1*(-Rd*xdot5 - Rd*xdot7 + Rl*xdot5 + Rl*xdot6 + 2*x1*xdot9 + x10*xdot2 - x11*xdot3 + x2*xdot10 - x3*xdot11 + 2*x9*xdot1) + K2*(-Rd*xdot5 - Rd*xdot7 - Rl*xdot5 + Rl*xdot6 + x1*xdot9 + 2*x10*xdot2 - x12*xdot4 + 2*x2*xdot10 - x4*xdot12 + x9*xdot1) + K3*(x1*xdot9 + x10*xdot2 + x11*xdot3 + x12*xdot4 + x2*xdot10 + x3*xdot11 + x4*xdot12 + x9*xdot1)]
    #~ ])

#~ F_eq = Fc

W1, W2, W3 = sp.symbols("W1, W2, W3", commutative=True)
# W1 = 2*L*(2*L*Rd + L*Rl + 2*Ll*Rd + Ll*Rl)
# W2 = (2*L**2*Rd + L**2*Rl + 2*L*Ll*Rd + 2*L*Ll*Rl + Ld*Ll*Rl)
# W3 = (2*L**2*Rd + 3*L**2*Rl + 2*L*Ld*Rl + 2*L*Ll*Rd + 2*L*Ll*Rl + Ld*Ll*Rl)
Fcalt = sp.Matrix([
[                                                                                                                                                                                                                                                                                                             L*xdot5 - L*xdot7 + Ll*xdot5 - Ll*xdot6 + Rl*x5 + Rl*x6 + Ul + x1*x9 - x11*x3],
        [                                                                                                                                                                                                                                                                                                            L*xdot6 - L*xdot8 - Ll*xdot5 + Ll*xdot6 - Rl*x5 + Rl*x6 - Ul + x10*x2 - x12*x4],
        [                                                                                                                                                                                                                                                                                                                                                                           C*xdot1 - x5*x9],
        [                                                                                                                                                                                                                                                                                                                                                                          C*xdot2 - x10*x6],
        [                                                                                                                                                                                                                                                                                                                                                                          C*xdot3 - x11*x7],
        [                                                                                                                                                                                                                                                                                                                                                                          C*xdot4 - x12*x8],
        [                                                                                                                                                                                                                                                                                                                                                      2*Rd*xdot5 + 2*Rd*xdot7 + 2*Rl*xdot6],
        [                                                                                                                                                                                                                                                                                                                                                             xdot5 - xdot6 + xdot7 - xdot8],
        [W1*(Rd*xdot5 + Rd*xdot7 + x11*xdot3 + x12*xdot4 + x3*xdot11 + x4*xdot12) + W2*(Rl*xdot5 + Rl*xdot6 + x1*xdot9 - x11*xdot3 - x3*xdot11 + x9*xdot1) + W3*(-Rl*xdot5 + Rl*xdot6 + x10*xdot2 - x12*xdot4 + x2*xdot10 - x4*xdot12)]])
F_eq = Fcalt

P10 = F_eq.jacobian(xxd)
#~ P10 = sp.Matrix([
        #~ [    0,      0,               0,               0,                L + Ll,           -Ll,    -L,  0,     0,     0,             0,             0],
        #~ [    0,      0,               0,               0,                   -Ll,        L + Ll,     0, -L,     0,     0,             0,             0],
        #~ [    C,      0,               0,               0,                     0,             0,     0,  0,     0,     0,             0,             0],
        #~ [    0,      C,               0,               0,                     0,             0,     0,  0,     0,     0,             0,             0],
        #~ [    0,      0,               C,               0,                     0,             0,     0,  0,     0,     0,             0,             0],
        #~ [    0,      0,               0,               C,                     0,             0,     0,  0,     0,     0,             0,             0],
        #~ [    0,      0,               0,               0,                  2*Rd,          2*Rl,  2*Rd,  0,     0,     0,             0,             0],
        #~ [    0,      0,               0,               0,                     1,            -1,     1, -1,     0,     0,             0,             0],
        #~ [W2*x9, W3*x10, W1*x11 - W2*x11, W1*x12 - W3*x12, Rd*W1 + Rl*W2 - Rl*W3, Rl*W2 + Rl*W3, Rd*W1,  0, W2*x1, W3*x2, W1*x3 - W2*x3, W1*x4 - W3*x4]])

P00 = F_eq.jacobian(xx)
#~ P00 = sp.Matrix([
        #~ [      x9,         0,                  -x11,                     0,  Rl,   Rl,    0,    0,       x1,        0,                 -x3,                   0],
        #~ [       0,       x10,                     0,                  -x12, -Rl,   Rl,    0,    0,        0,       x2,                   0,                 -x4],
        #~ [       0,         0,                     0,                     0, -x9,    0,    0,    0,      -x5,        0,                   0,                   0],
        #~ [       0,         0,                     0,                     0,   0, -x10,    0,    0,        0,      -x6,                   0,                   0],
        #~ [       0,         0,                     0,                     0,   0,    0, -x11,    0,        0,        0,                 -x7,                   0],
        #~ [       0,         0,                     0,                     0,   0,    0,    0, -x12,        0,        0,                   0,                 -x8],
        #~ [       0,         0,                     0,                     0,   0,    0,    0,    0,        0,        0,                   0,                   0],
        #~ [       0,         0,                     0,                     0,   0,    0,    0,    0,        0,        0,                   0,                   0],
        #~ [W2*xdot9, W3*xdot10, W1*xdot11 - W2*xdot11, W1*xdot12 - W3*xdot12,   0,    0,    0,    0, W2*xdot1, W3*xdot2, W1*xdot3 - W2*xdot3, W1*xdot4 - W3*xdot4]])

P10_roc = st.nullspaceMatrix(P10)
#~ P10_roc = sp.Matrix([
        #~ [     0,             0,             0],
        #~ [     0,             0,             0],
        #~ [     0,             0,             0],
        #~ [     0,             0,             0],
        #~ [     0,             0,             0],
        #~ [     0,             0,             0],
        #~ [     0,             0,             0],
        #~ [     0,             0,             0],
        #~ [-W3*x2, x3*(-W1 + W2), x4*(-W1 + W3)],
        #~ [ W2*x1,             0,             0],
        #~ [     0,         W2*x1,             0],
        #~ [     0,             0,         W2*x1]])

# statt P10_rpinv nehme ich jetzt S mit einheitsvektoren sodass
# (S,P10_roc) regulär ist:



#~ from IPython import embed as IPS
#~ IPS()
    
if 0:    
        
    # algebraische Substitutionsgleichungen mit Ableitungen:
    x6subs = -Rd*(x5 + x7)/Rl
    xdot6subs = -Rd*(xdot5 + xdot7)/Rl
    x8subs =  (x5 + x7)*(Rd/Rl + 1)
    xdot8subs = (xdot5 + xdot7)*(Rd/Rl + 1)

    F_eq_new = F_eq.subs(x8,x8subs).subs(xdot8,xdot8subs).subs(x6,x6subs).subs(xdot6,xdot6subs)

    P10 = F_eq_new.jacobian(vec_xdot)

