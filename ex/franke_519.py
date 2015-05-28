# -*- coding: utf-8 -*-
import sympy as sp

m0, m1, m2 = sp.symbols("m0, m1, m2")
l0, l1, l2 = sp.symbols("l0, l1, l2")
J0, J1, J2 = sp.symbols("J0, J1, J2")
L1, g = sp.symbols("L1, g")

x1, x2, x3, x4, x5, x6 = sp.symbols("x1, x2, x3, x4, x5, x6")
vec_x = sp.Matrix([x1, x2, x3, x4, x5, x6])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4, xdot5, xdot6 = vec_xdot




# Abkuerzungen:
a2 = J0 + (m0+m1+m2)*l0**2 + (m1*l1 + m2*L1)*l0*sp.cos(x2) + m2*l2*l0*sp.cos(x3)
b2 = (m1*l1 + m2*L1)*l0*sp.cos(x2) + J1 + m1*l1**2 + m2*L1**2 + m2*l2*L1*sp.cos(-x3+x2)
c2 = m2*l2*l0*sp.cos(x3) + m2*l2*L1*sp.cos(-x3+x2) + J2 + m2*l2**2
b1 = -2*l0*(m1*l1 + m2*L1)*sp.sin(x2)*x5 + 2*m2*l1*L1*sp.sin(x3-x2)*x5
c1 = -2*m2*l2*l0*sp.sin(x3)*x6 + 2*m2*l2*L1*sp.sin(x3-x2)*x6
b0 = (m1*l1 + m2*L1)*( (-g-l0*x5**2)*sp.cos(x2) - l0*sp.sin(x2)*(xdot4 + xdot5) ) + \
        m2*l2*L1*( sp.sin(x3-x2)*(xdot6+xdot5) + sp.cos(x3-x2)*(x6**2-x5**2) )
c0 = m2*l2*( (-g-l0*x6**2)*sp.cos(x3) - l0*sp.sin(x3)*(xdot4+xdot6) ) + \
        m2*l2*L1*( -sp.sin(x3-x2)*(xdot6 + xdot5) - sp.cos(x3-x2)*(x6**2 - x5**2) )




#diffvec_x = sp.Matrix([x1, x2, x3, x4, x5, x6, xdot1, xdot2, xdot3, xdot4, xdot5, xdot6 ])




P1i = sp.Matrix([ [1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,a2,b2,c2] ])
P0i = sp.Matrix([ [0,0,0,-1,0,0],[0,0,0,0,-1,0],[0,0,0,0,0,-1],[0,b0,c0,0,b1,c1] ])
