# -*- coding: utf-8 -*-
import sympy as sp
import rst_symbtools.symb_tools as st

Nz = 5  # Anzahl der Zustandskomponenten
Nin = 2  # Anzahl der Eingänge

Neq = Nz - Nin  # Anzahl der Systemgleichungen


eq_coeffsA = []
eq_coeffsB = []
#eq_coeffsC = []

all_coeffs = []
# iterieren über die Gleichungen
for i in xrange(Neq):
    coeffsA = []
    coeffsB = []
    #coeffsC = []
    coeffs = []
    # iterieren über die einzelnen Freiheitsgrade
    for j in xrange(Nz):
        # drei Ableitungsordnungen
        A = sp.Symbol("A%i_%i" % (i + 1, j + 1))
        B = sp.Symbol("B%i_%i" % (i + 1, j + 1))
        #C = sp.Symbol("C%i_%i" % (i + 1, j + 1))
        
        coeffsA.append(A)
        coeffsB.append(B)
        #coeffsC.append(C)
        
        coeffs.extend([A, B])
     
     # jede Gleichung bekommt eine separate Liste
    eq_coeffsA.append(coeffsA)
    eq_coeffsB.append(coeffsB)
    #eq_coeffsC.append(coeffsC)
     
    all_coeffs.extend(coeffs)
        

diff_symbols = sp.Matrix(all_coeffs)

#x1, x2, x3, x4, x5, x6 = sp.symbols("x1, x2, x3, x4, x5, x6")


xx = vec_x = sp.Matrix( sp.symbols("x1:%i" % (Nz + 1)) )

vec_xdot = st.perform_time_derivative(vec_x, vec_x)

AA = sp.Matrix(eq_coeffsA)
BB = sp.Matrix(eq_coeffsB)
BB = sp.zeros(Neq, Nz)
BB[:, :Neq] = sp.eye(Neq)
BB[:, Neq:] = sp.Matrix([[1, 0], [1, 0], [0, 1]])


eq_sys = AA*xx + BB*vec_xdot

F_eq = eq_sys




