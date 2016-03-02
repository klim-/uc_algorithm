# -*- coding: utf-8 -*-
import sympy as sp
import symbtools as st

Ndf = 4  # Anzahl der Freiheitsgrade
Nin = 3  # Anzahl der Eingänge

Neq = Ndf - Nin  # Anzahl der mechanischen Gleichungen (ohne definitorische)


eq_coeffsA = []
eq_coeffsB = []
eq_coeffsC = []

all_coeffs = []
# iterieren über die Gleichungen
for i in xrange(Neq):
    coeffsA = []
    coeffsB = []
    coeffsC = []
    coeffs = []
    # iterieren über die einzelnen Freiheitsgrade
    for j in xrange(Ndf):
        # drei Ableitungsordnungen
        A = sp.Symbol("A%i_%i" % (i + 1, j + 1))
        B = sp.Symbol("B%i_%i" % (i + 1, j + 1))
        C = sp.Symbol("C%i_%i" % (i + 1, j + 1))
        
        coeffsA.append(A)
        coeffsB.append(B)
        coeffsC.append(C)
        
        coeffs.extend([A, B, C])
     
     # jede Gleichung bekommt eine separate Liste
    eq_coeffsA.append(coeffsA)
    eq_coeffsB.append(coeffsB)
    eq_coeffsC.append(coeffsC)
     
    all_coeffs.extend(coeffs)
        

diff_symbols = sp.Matrix(all_coeffs)

#x1, x2, x3, x4, x5, x6 = sp.symbols("x1, x2, x3, x4, x5, x6")


vec_x = sp.Matrix( sp.symbols("x1:%i" % (Ndf*2 + 1)) )

theta = vec_x[:Ndf, :]
mu = vec_x[Ndf:, :]

vec_xdot = st.time_deriv(vec_x, vec_x)
thetadot = st.time_deriv(theta, vec_x)
mudot = st.time_deriv(mu, vec_x)


# definitorische Gleichungen
eq_defin = thetadot - mu

AA = sp.Matrix(eq_coeffsA)
BB = sp.Matrix(eq_coeffsB)
CC = sp.Matrix(eq_coeffsC)


#eq_mech = AA*theta + BB*mu + CC*mudot
eq_mech = AA*theta + BB*thetadot + CC*mudot

F_eq = st.row_stack(eq_defin, eq_mech)

#F_eq = sp.Matrix([
        #[ xdot1 - x4 ],
        #[ xdot2 - x5 ],
        #[ xdot3 - x6 ],
        #[ a2*xdot4 + b2*xdot5 + c2*xdot6 + b1*x5 + c1*x6 + b0*x2 + c0*x3 ]])
        #~ [ a2*xdot4 + b2*xdot5 + c2*xdot6 + b1*x5 + c1*x6 + b0*x2 + c0*x3 + a0*x1 + a1*x4]])
