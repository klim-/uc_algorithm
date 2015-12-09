# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
from sympy import sin, cos
import symb_tools as st

"Beispiel aus 2003_Rouchon_Murray_Martin_TechReport_Flat_systems_equivalence_trajectory_generation.pdf"



a = sp.symbols("a")

# x = theta, phi, psi omega
x1, x2, x3, x4 = xx =  st.symb_vector("x1:5")
vec_x = xx

vec_xdot = xx_dot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4 = xx_dot

eq1 = xdot4 - a*(xdot1 - xdot3*sin(x2))*( (xdot2 - x4*sin(x1))/cos(x1))
eq2 = xdot3 - (x4 - xdot2*sin(x1))/(cos(x2)*cos(x1))
F_eq_orig = sp.Matrix([ eq1, eq2])


# Algorithmus läuft nicht durch (assertion error)

P0_orig = F_eq_orig.jacobian(xx)
P1_orig = F_eq_orig.jacobian(xx_dot)



P0, replm_A, tdsA  = st.introduce_abreviations(P0_orig, prefix="A", time_dep_symbs=xx)
P1, replm_B, tdsB  = st.introduce_abreviations(P1_orig, prefix="B", time_dep_symbs=xx)


subs_tuples = replm_A + replm_B
diff_symbols = st.row_stack(tdsA, tdsB)

F_eq = P0*xx + P1*xx_dot

from IPython import embed as IPS
IPS()