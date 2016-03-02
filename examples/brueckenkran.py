# -*- coding: utf-8 -*-
import sympy as sp
import symbtools as st

# Differentialgleichung 2. Ordnung:
#~ sp.Matrix([
        #~ [m2*q2*(g*sin(p1) + pddot1*q2 + 2*pdot1*qdot2 + qddot1*cos(p1))]
    #~ ])

n = 3  # Anzahl der Freiheitsgrade
nq = 2  # Anzahl der Eingänge

Neq = n - nq  # Anzahl der mechanischen Gleichungen (ohne definitorische)

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
    for j in xrange(n):
        # drei Ableitungsordnungen
        A = sp.Symbol("A%i%i" % (i + 1, j + 1))
        B = sp.Symbol("B%i%i" % (i + 1, j + 1))
        C = sp.Symbol("C%i%i" % (i + 1, j + 1))
        
        coeffsA.append(A)
        coeffsB.append(B)
        coeffsC.append(C)
        
        coeffs.extend([A, B, C])
     
     # jede Gleichung bekommt eine separate Liste
    eq_coeffsA.append(coeffsA)
    eq_coeffsB.append(coeffsB)
    eq_coeffsC.append(coeffsC)
     
    all_coeffs.extend(coeffs)


vec_x = sp.Matrix( sp.symbols("x1:%i" % (n*2 + 1)) )

theta = vec_x[:n, :]
mu = vec_x[n:, :]

vec_xdot = st.time_deriv(vec_x, vec_x)
thetadot = st.time_deriv(theta, vec_x)
mudot = st.time_deriv(mu, vec_x)



AA = sp.Matrix(eq_coeffsA) # P0
BB = sp.Matrix(eq_coeffsB) # P1
CC = sp.Matrix(eq_coeffsC) # P2


# Informationen über das konkrete System:
diff_symbols = [AA[0,0], AA[0, 2], BB[0, 0], BB[0, 2], CC[0, 0], CC[0, 1]]

diff_symbols = sp.Matrix(sorted(diff_symbols, key=str))

AA[0,1]*=0

BB[0,1]*=0

CC[0,2]*=0

# definitorische Gleichungen
eq_defin = thetadot - mu

#~ eq_mech = AA*theta + BB*mu + CC*mudot
eq_mech = AA*theta + BB*thetadot + CC*mudot

# Container um zusätzliche Information über das Beispiel zu speichern
data = st.Container()
data.P0 = AA
data.P1 = BB
data.P2 = CC
data.eq_mech = eq_mech
data.time_dep_symbols = diff_symbols


F_eq_orig = F_eq = st.row_stack(eq_defin, eq_mech)
