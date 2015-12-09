# -*- coding: utf-8 -*-
import sympy as sp
import rst_symbtools.symb_tools as st
import itertools

from sympy import sin, cos

Ndf = 2  # Anzahl der Freiheitsgrade
Nin = 1  # Anzahl der Eingänge

Neq = Ndf - Nin  # Anzahl der mechanischen Gleichungen (ohne definitorische)


vec_x = sp.Matrix( sp.symbols("x1:%i" % (Ndf*2 + 1)) )
vec_xdot = st.perform_time_derivative(vec_x, vec_x)
st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)

params = sp.symbols('m1, m2, s2, g')
st.make_global(params, 1)




F_eq = sp.Matrix([[x3 - xdot1], [x4 - xdot2], [-g*sin(x1)/s2 - xdot3 - xdot4*cos(x1)/s2]])
# Diese Gleichung wird durch die Einführung allgemeiner Zeitfunktionen
# ("zeitabhägiger Symbole") vereinfacht.

P1 = F_eq.jacobian(vec_xdot)
P0 = F_eq.jacobian(vec_x)

P1_new = P1*1
P0_new = P0*1

Nc, Nr = P1.shape


gen0 = sp.numbered_symbols('A')
gen1 = sp.numbered_symbols('B')

all_symbols = []
diff_symbols = []
subs_tuples = []

for i,j in itertools.product(range(Nc), range(Nr)):
    elt0 = P0_new[i, j]
    if not st.is_number(elt0):
        symb = gen0.next()
        all_symbols.append(symb)
        P0_new[i,j] = symb
        subs_tuples.append((symb, elt0))
        if st.depends_on_t(elt0, 't', vec_x):
            diff_symbols.append(symb)
        
    elt1 = P1_new[i, j]
    if not st.is_number(elt1):
        symb = gen1.next()
        all_symbols.append(symb)
        P1_new[i,j] = symb
        subs_tuples.append((symb, elt1))
        if st.depends_on_t(elt1, 't', vec_x):
            diff_symbols.append(symb)
    #print i, j, elt0, elt1
        

F_eq = P1_new*vec_xdot + P0_new*vec_x

diff_symbols = sp.Matrix(diff_symbols)

from IPython import embed as IPS
IPS()