# -*- coding: utf-8 -*-
# enable true divison
#from __future__ import division
import sympy as sp
import symbtools as st

I1, I2 = sp.symbols("I1, I2")
m1, m2 = sp.symbols("m1, m2")
l, r = sp.symbols("l, r")
g = sp.symbols("g")


x1, x2, x3, x4 = sp.symbols("x1, x2, x3, x4")
vec_x = sp.Matrix([x1, x2, x3, x4])

vec_xdot = st.perform_time_derivative(vec_x, vec_x)
xdot1, xdot2, xdot3, xdot4 = vec_xdot

F_eq = sp.Matrix([
        [xdot1 - x2],
        [xdot2 + (((m1*r+m2*l)*g)/(I1+m2*l**2))*sp.sin(x1) + ((I2)/(I1+m2*l**2))*xdot4],
        [xdot3-x4]])


#P1i = sp.Matrix([ [1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0] ])
#P0i = sp.Matrix([ [0,0,-x4,-x3,0],[0,0,0,-1,0],[0,0,0,0,-1] ])
