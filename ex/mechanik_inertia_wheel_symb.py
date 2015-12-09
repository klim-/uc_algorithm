# -*- coding: utf-8 -*-
import sympy as sp
import rst_symbtools.symb_tools as st
import itertools

from sympy import sin, cos

"""
Acrobot Modell auf Basis der implizieten Systemgleichungen 2. Ordnung
"""



Ndf = 2  # Anzahl der Freiheitsgrade
Nin = 1  # Anzahl der Eingänge

Neq = Ndf - Nin  # Anzahl der mechanischen Gleichungen (ohne definitorische)


vec_x = sp.Matrix( sp.symbols("x1:%i" % (Ndf*2 + 1)) )
vec_xdot = st.perform_time_derivative(vec_x, vec_x)
st.make_global(vec_x, 1)
st.make_global(vec_xdot, 1)

P0 = sp.sympify("""
Matrix([
[ 0, 0, 1, 0],
[ 0, 0, 0, 1],
[A0, 0, 0, 0]])""")

P1 = sp.sympify("""
Matrix([
[-1,  0,  0,  0],
[ 0, -1,  0,  0],
[ 0,  0, C0, C1]])""")

st.make_global(P0.s, 1)
st.make_global(P1.s, 1)

F_eq = P0*vec_x + P1*vec_xdot
diff_symbols = sp.Matrix(
[
[A0]])


from IPython import embed as IPS
IPS()