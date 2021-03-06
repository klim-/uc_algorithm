# -*- coding: utf-8 -*-
import sympy as sp
import symbtools as st

n = 4  # Anzahl der Freiheitsgrade
nq = 2  # Anzahl der Eingänge

Neq = n - nq  # Anzahl der mechanischen Gleichungen (ohne definitorische)


eq_coeffsA = []
eq_coeffsB = []
eq_coeffsC = []

all_coeffs = []
# iterieren über die Gleichungen
for i in range(Neq):
    coeffsA = []
    coeffsB = []
    coeffsC = []
    coeffs = []
    # iterieren über die einzelnen Freiheitsgrade
    for j in range(n):
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
        



#x1, x2, x3, x4, x5, x6 = sp.symbols("x1, x2, x3, x4, x5, x6")


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
diff_symbols = list(AA[:,:2]) + list(AA[:, -1]) + list(BB[:, 2])

diff_symbols = sp.Matrix(sorted(diff_symbols, key=str))

AA[:,2]*=0

BB[:,:2]*=0
BB[:,3]*=0

CC[:,2]*=0

CC[:,3]=CC[:,0]
CC[0,1]=CC[1,0]
CC[1,1]=CC[1,0]

# definitorische Gleichungen
eq_defin = thetadot - mu

#eq_mech = AA*theta + BB*mu + CC*mudot
eq_mech = AA*theta + BB*thetadot + CC*mudot

# Container um zusätzliche Information über das Beispiel zu speichern
data = st.Container()
data.P0 = AA
data.P1 = BB
data.P2 = CC
data.eq_mech = eq_mech
data.time_dep_symbols = diff_symbols


F_eq_orig = F_eq = st.row_stack(eq_defin, eq_mech)


sys_name = "mechanik_ph_RTtt_seriell"



#~ from ipHelp import IPS
#~ IPS()


#F_eq = sp.Matrix([
        #[ xdot1 - x4 ],
        #[ xdot2 - x5 ],
        #[ xdot3 - x6 ],
        #[ a2*xdot4 + b2*xdot5 + c2*xdot6 + b1*x5 + c1*x6 + b0*x2 + c0*x3 ]])
        #~ [ a2*xdot4 + b2*xdot5 + c2*xdot6 + b1*x5 + c1*x6 + b0*x2 + c0*x3 + a0*x1 + a1*x4]])
