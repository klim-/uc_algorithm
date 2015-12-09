# -*- coding: utf-8 -*-
import sympy as sp
import rst_symbtools.symb_tools as st

Ndf = 4  # Anzahl der Freiheitsgrade
Nin = 2  # Anzahl der Eingänge

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
        



vec_x = sp.Matrix( sp.symbols("x1:%i" % (Ndf*2 + 1)) )

theta = vec_x[:Ndf, :]
mu = vec_x[Ndf:, :]

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
thetadot = st.perform_time_derivative(theta, vec_x)
mudot = st.perform_time_derivative(mu, vec_x)


# definitorische Gleichungen
eq_defin = thetadot - mu

st.make_global(all_coeffs, 1)



# Weil die allgemeine Struktur so kompliziert ist, dass der Algorithmus zu lange
# braucht, werden nur die Terme (mit Platzhaltern) berücksichtigt,
# die tatsächlich auftreten
AA = sp.Matrix([[A1_1, 0, 0, 0], [0, A2_2, 0, 0]])
BB = sp.zeros(Neq, Ndf)
CC = sp.Matrix([[C1_1, 0, C1_3, C1_4], [0, C2_2, C2_3, 0]])


# C1_1, C2_2 sind constant
diff_symbols = sp.Matrix([A1_1, C1_3, C1_4])


#eq_mech = AA*theta + BB*mu + CC*mudot
eq_mech = AA*theta + BB*thetadot + CC*mudot



F_eq = st.row_stack(eq_defin, eq_mech)


from IPython import embed as IPS
IPS()
