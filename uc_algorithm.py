# -*- coding: utf-8 -*-

#    Copyright (C) 2015-2016
#    by Klemens Fritzsche, 2e63a67d46@leckstrom.de
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Implementation of an algorithm for determining differential flat ouputs
for nonlinear dynamic systems (see Franke/Röbenack: On the Computation
of Flat Outputs for Nonlinear Control Systems. European Control
Conference [ECC], 2013).
"""

# enable true divison
from __future__ import division
import sys
import importlib

from IPython import embed as IPS

import numpy as np
import sympy as sp
import symbtools as st
import symbtools.noncommutativetools as nct

import core.algebra as al
import core.matrix_container as mc
import core.system_container as sc
import core.integrability as ic
import core.transformations as tr
import util.print_candy as pc


try:
    # import example passed by command line argument
    path = sys.argv[1]
    pseudo_path = path.replace("/",".").replace(".py","")
    example = importlib.import_module(pseudo_path)
except ImportError:
    raise ImportError('The argument ' + str(sys.argv[1]) + ' does not seem to point to a valid example.')

# vector for optional time dependent variables, specified in examples file
if hasattr(example, 'diff_symbols'):
    diff_symbols = example.diff_symbols
else:
    diff_symbols = sp.Matrix([])

mode = "auto" # "manual" or "auto"

calc_G = False # calculate G matrix

# raw_input() does not seem to handle empty strings, so this is a
# workaround to escape from the input
auto = None

myStack = sc.SystemStack(diff_symbols)
myStack.vec_x = example.vec_x
myStack.calc_G = calc_G
print "x ="; pc.print_nicely(myStack.vec_x)
print "\n\n xdot ="; pc.print_nicely(myStack.vec_xdot)

def end_condition(Bi):
    """ The algorithm ends if Bi has full row rank.
    """
    n, p = Bi.shape
    return True if (n == st.rnd_number_rank(Bi)) else False

def is_special_case(Bi):
    """ Checks for special case
    """
    n, p = Bi.shape
    return True if (st.rnd_number_rank(Bi) < p) else False

def roc_hint(matrix_string, i,  matrix):
    error_string = "There must have been a mistake. Try again.\n"
    pc.print_matrix(matrix_string, i, "", matrix)
    while True:
        try:
            roc = eval(raw_input("Please enter " + str(matrix_string) + str(i) + "_roc or \"auto\":\n"))
            if roc==auto:
                roc = al.right_ortho_complement(matrix)
            try:
                if al.is_zero_matrix(matrix*roc):
                    pc.print_matrix(matrix_string, i, "_roc", roc)
                    return roc
                else:
                    print error_string
            except:
                print error_string
        except Exception as exc:
            print exc
            print error_string

def loc_hint(matrix_string, i,  matrix):
    error_string = "There must have been a mistake. Try again.\n"
    pc.print_matrix(matrix_string, i, "", matrix)
    while True:
        try:
            loc = eval(raw_input("Please enter " + str(matrix_string) + str(i) + "_loc or \"auto\":\n"))
            if loc==auto:
                loc = al.left_ortho_complement(matrix)
            try:
                if al.is_zero_matrix(loc*matrix):
                    pc.print_matrix(matrix_string, i, "_loc", loc)
                    return loc
                else:
                    print error_string
            except:
                print error_string
        except Exception as exc:
            print exc
            print error_string

def rpinv_hint(matrix_string, i,  matrix):
    # TODO: check if this can somehow move to the algebra module
    #       (problem: different namespace, so the algebra module won't
    #                 know about symbols and statevector)
    error_string = "There must have been a mistake. Try again.\n"
    while True:
        try:
            rpinv = eval(raw_input("Please enter " + str(matrix_string) + str(i) + "_rpinv or \"auto\":\n"))
            if rpinv==auto:
                rpinv = al.right_pseudo_inverse(matrix)
            try:
                if al.is_unit_matrix(matrix*rpinv):
                    pc.print_matrix(matrix_string, i, "_rpinv", rpinv)
                    return rpinv
                else:
                    print error_string
            except:
                print error_string
        except Exception as exc:
            print exc
            print error_string

def lpinv_hint(matrix_string, i,  matrix):
    # TODO: check if this can somehow move to the algebra module
    #       (problem: different namespace, so the algebra module won't
    #                 know about symbols and statevector)
    error_string = "There must have been a mistake. Try again.\n"
    while True:
        try:
            lpinv = eval(raw_input("Please enter " + str(matrix_string) + str(i) + "_lpinv or \"auto\":\n"))
            if lpinv==auto:
                lpinv = al.left_pseudo_inverse(matrix)
            try:
                if al.is_unit_matrix(lpinv*matrix):
                    pc.print_matrix(matrix_string, i, "_lpinv", lpinv)
                    return lpinv
                else:
                    print error_string
            except:
                print error_string
        except Exception as exc:
            print exc
            print error_string

def reduction(iter_stack):
    P1i = iter_stack.P1
    P0i = iter_stack.P0
    
    P1i_roc = roc_hint("P1", iter_stack.i, P1i) if mode=="manual" else al.right_ortho_complement(P1i)
    P1i_rpinv = rpinv_hint("P1", iter_stack.i, P1i) if mode=="manual" else al.right_pseudo_inverse(P1i)

    P1i_dot = st.time_deriv(P1i, myStack.diffvec_x)

    Ai = al.custom_simplify( (P0i - P1i_dot)*P1i_rpinv )
    Bi = al.custom_simplify( (P0i - P1i_dot)*P1i_roc )

    if myStack.calc_G:
        Bi_lpinv = al.left_pseudo_inverse(Bi)
    else:
        Bi_lpinv = None

    # store
    iter_stack.store_reduction_matrices( Ai, Bi, Bi_lpinv, P1i_roc, P1i_rpinv, P1i_dot )

    return Bi

def fourseven(iter_stack):
    # load matrices
    Ai = iter_stack.A
    Bi = iter_stack.B
    P1i_roc = iter_stack.P1_roc
    P1i_rpinv = iter_stack.P1_rpinv
    
    K2 = roc_hint("B", iter_stack.i, Bi) if mode=="manual" else al.right_ortho_complement(Bi)

    K1 = al.regular_completion(K2)

    K = st.concat_cols(K1, K2)

    assert al.is_regular_matrix(K), "K is not a regular matrix."

    Bi_tilde = al.custom_simplify(Bi*K1) # unit matrix

    if myStack.calc_G:
        Bi_tilde_lpinv = al.left_pseudo_inverse(Bi_tilde)
    else:
        Bi_tilde_lpinv = None

    P1i_tilde_roc = al.custom_simplify( P1i_roc*K1 )

    Zi = al.custom_simplify( P1i_roc*K2 )

    Zi_lpinv = al.Zi_left_pinv_with_restrictions(P1i_rpinv, P1i_tilde_roc, Zi)

    assert al.is_unit_matrix( Zi_lpinv*Zi ), "Zi_lpinv seems to be wrong."

    # store
    iter_stack.store_special_case_matrices( Zi, Zi_lpinv, Bi_tilde, Bi_tilde_lpinv, P1i_tilde_roc )

    return Bi_tilde

def elimination(iter_stack):
    # load matrices
    Ai = iter_stack.A
    Bi = iter_stack.B
    
    Bi_loc = loc_hint("B", iter_stack.i, Bi) if mode=="manual" else al.left_ortho_complement(Bi)
        
    # store
    iter_stack.store_elimination_matrices(Bi_loc)

    P1i_new = Bi_loc
    P0i_new = Bi_loc*Ai

    return P1i_new, P0i_new

def tangent_system():
    # exterior derivative of F_eq:
    try:
        P1i # in case the system is given by the matrices P1i and P0i
    except NameError:
        print "\n\n0 = F(x,xdot) ="; pc.print_nicely(example.F_eq)

        P1i = example.F_eq.jacobian(myStack.vec_xdot)
        P0i = example.F_eq.jacobian(myStack.vec_x)

    print "\n\n"
    return P1i, P0i
    

########################################################################

def main():
    P1i, P0i = tangent_system()

    global mode
    mode = raw_input("Enter \"auto\" for automatic mode, or \"manual\" for manual mode:\n")
    #####################################################################

    i = 0

    while 1:
        # new iteration_stack
        myIteration = mc.IterationStack(i, P1i, P0i)

        assert al.has_full_row_rank(P1i), "P10 does not have full row rank. There \
                                        must be algebraic equations."

        Bi = reduction(myIteration)

        assert not al.is_zero_matrix(Bi), "System ist not flat!"

        if is_special_case(Bi):
            # special case
           Bi = fourseven(myIteration)

        if end_condition(Bi):
            # flag
            myIteration.last_iter_step = True

            # add to system stack + print
            myStack.add_iteration(myIteration)
            break

        P1i, P0i = elimination(myIteration)

        if mode=="manual": pc.print_line()

        # add to system stack + print iteration
        myStack.add_iteration(myIteration)

        i += 1

    # create transformation and calculate Q and G(d/dt)
    myStack.transformation = tr.Transformation(myStack)


    # store results of the algorithm in pickled data container
    #~ if hasattr(example, 'data'):
        #~ data = example.data
    data = st.Container()
    data.F_eq = nct.make_all_symbols_noncommutative(example.F_eq, "")[0]
    data.vec_x = nct.make_all_symbols_noncommutative(example.vec_x, "")[0]
    data.vec_xdot = nct.make_all_symbols_noncommutative(example.vec_xdot, "")[0]
    data.diff_symbols = nct.make_all_symbols_noncommutative(example.diff_symbols, "")[0]

    # add data to be stored here:
    # make symbols in P10 and P00 noncommutative
    P1 = nct.make_all_symbols_noncommutative(myStack.transformation.P10, "")[0]
    P0 = nct.make_all_symbols_noncommutative(myStack.transformation.P00, "")[0]
    s  = sp.Symbol('s', commutative=False)
    data.P = P1*s+P0
    data.P1 = P1
    data.P0 = P0
    data.Q = nct.make_all_symbols_noncommutative(myStack.transformation.Q, "")[0]

    fname = path.replace(".py",".pcl")
    st.pickle_full_dump(data, fname)

    # for testing
    try:
        import pycartan as ct
        global w
        myIntegrabilityCheck = ic.IntegrabilityCheck(myStack)
        w = myStack.transformation.w
    except Exception as exc:
        print exc

    print "Data saved to ", fname, "\n"
    pc.print_line()

if __name__ == '__main__':
    main()

