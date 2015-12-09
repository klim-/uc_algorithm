# -*- coding: utf-8 -*-
# enable true divison
from __future__ import division
import sympy as sp
from sympy import sin, cos
import symb_tools as st


# x = theta, phi, psi omega
x1, x2, x3, x4 = xx =  st.symb_vector("x1:5")

vec_xdot = xx_dot = st.perform_time_derivative(xx, xx)
xdot1, xdot2, xdot3, xdot4 = xx_dot

vec_x = xx


nx = len(xx)
nu = 2
nF = nx - nu


P0 = st.symbMatrix(nF, nx, 'A')
P1 = st.symbMatrix(nF, nx, 'B')


diff_symbols = sp.Matrix(list(P0) + list(P1))

F_eq_orig = F_eq = P0*xx + P1*xx_dot

from IPython import embed as IPS
IPS()



