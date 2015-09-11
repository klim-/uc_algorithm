# -*- coding: utf-8 -*-

"""
Erste Implementierung des Algorithmus zur Konstruktion flacher Ausgänge für
nichtlineare Systeme von Dr. Matthias Franke, TU Dresden
"""

# enable true divison
from __future__ import division
import sys

from IPython import embed as IPS

import numpy as np
import sympy as sp
import symb_tools as st
import diffgeopy as ct


from algebra import *
from integrability import *
from matrix_container import *
from print_candy import *
from system_container import *
from transformations import *


from examples import *


mode = "auto" # "manual" or "auto"

# raw_input() does not seem to handle empty strings, so this is a
# workaround to escape from the input
auto = None

myStack = SystemStack(diff_symbols)
myStack.vec_x = vec_x
print "x ="; print_nicely(myStack.vec_x)
print "\n\n xdot ="; print_nicely(myStack.vec_xdot)

def end_condition(Bi):
    """ The algorithm ends if Bi has full row rank.
    """
    n, p = Bi.shape
    return True if (n == srank(Bi)) else False

def outlier(Bi):
    """ sonderfall 4.7
        tritt auf wenn Bi linear abhängige spalten aufweisst bzw.
        keinen vollen spaltenrang hat
    """
    n, p = Bi.shape
    return True if (srank(Bi) < p) else False

def roc_hint(matrix_string, i,  matrix):
    error_string = "There must have been a mistake. Try again.\n"
    print_matrix(matrix_string, i, "", matrix)
    while True:
        try:
            roc = eval(raw_input("Please enter " + str(matrix_string) + str(i) + "_roc or \"auto\":\n"))
            if roc==auto:
                roc = right_ortho_complement(matrix)
            try:
                if is_zero_matrix(matrix*roc):
                    print_matrix(matrix_string, i, "_roc", roc)
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
    print_matrix(matrix_string, i, "", matrix)
    while True:
        try:
            loc = eval(raw_input("Please enter " + str(matrix_string) + str(i) + "_loc or \"auto\":\n"))
            if loc==auto:
                loc = left_ortho_complement(matrix)
            try:
                if is_zero_matrix(loc*matrix):
                    print_matrix(matrix_string, i, "_loc", loc)
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
                rpinv = right_pseudo_inverse(matrix)
            try:
                if is_unit_matrix(matrix*rpinv):
                    print_matrix(matrix_string, i, "_rpinv", rpinv)
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
                lpinv = left_pseudo_inverse(matrix)
            try:
                if is_unit_matrix(lpinv*matrix):
                    print_matrix(matrix_string, i, "_lpinv", lpinv)
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
    
    P1i_roc = roc_hint("P1", iter_stack.i, P1i) if mode=="manual" else right_ortho_complement(P1i)
    P1i_rpinv = rpinv_hint("P1", iter_stack.i, P1i) if mode=="manual" else right_pseudo_inverse(P1i)


    P1i_dot = st.perform_time_derivative(P1i, myStack.diffvec_x)

    Ai = custom_simplify( (P0i - P1i_dot)*P1i_rpinv )
    Bi = custom_simplify( (P0i - P1i_dot)*P1i_roc )

    # TODO: this is not very elegant!
    try:
        Bi_lpinv = left_pseudo_inverse(Bi)
    except:
        Bi_lpinv = None

    # store
    iter_stack.store_reduction_matrices( Ai, Bi, Bi_lpinv, P1i_roc, P1i_rpinv, P1i_dot )

    return Bi

def fourseven(iter_stack):#, Ai, Bi, P1i_roc, P1i_rpinv):
    # load matrices
    Ai = iter_stack.A
    Bi = iter_stack.B
    P1i_roc = iter_stack.P1_roc
    P1i_rpinv = iter_stack.P1_rpinv
    
    K2 = roc_hint("B", iter_stack.i, Bi) if mode=="manual" else right_ortho_complement(Bi)


    if has_full_row_rank(Bi):
        K1 = rpinv_hint("B", iter_stack.i, Bi) if mode=="manual" else right_pseudo_inverse(Bi)
    else:
        K1 = roc_hint("K2.T", iter_stack.i, K2.T) if mode=="manual" else right_ortho_complement(K2.T)

    K = st.concat_cols(K1, K2)

    assert is_regular_matrix(K), "K is not a regular matrix."

    Bi_tilde = custom_simplify(Bi*K1) # unit matrix
    Bi_tilde_lpinv = left_pseudo_inverse(Bi_tilde)

    P1i_tilde_roc = custom_simplify( P1i_roc*K1 )

    Zi = custom_simplify( P1i_roc*K2 )

    Zi_lpinv = Zi_left_pinv_with_restrictions(P1i_rpinv, P1i_tilde_roc, Zi)

    assert is_unit_matrix( Zi_lpinv*Zi ), "Zi_lpinv seems to be wrong."

    # store
    iter_stack.store_outlier_matrices( Zi, Zi_lpinv, Bi_tilde, Bi_tilde_lpinv, P1i_tilde_roc )

    return Bi_tilde

def elimination(iter_stack):
    # load matrices
    Ai = iter_stack.A
    Bi = iter_stack.B
    
    Bi_loc = loc_hint("B", iter_stack.i, Bi) if mode=="manual" else left_ortho_complement(Bi)
        
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
        print "\n\n0 = F(x,xdot) ="; print_nicely(F_eq)

        P1i = sp.Matrix([])
        P0i = sp.Matrix([])
        for i in xrange(len(myStack.vec_xdot)):
            vector1 = custom_simplify( F_eq.diff(myStack.vec_xdot[i]) )
            P1i = st.concat_cols(P1i, vector1)

            vector0 = custom_simplify( F_eq.diff(myStack.vec_x[i]) )
            P0i = st.concat_cols(P0i, vector0)
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
        myIteration = IterationStack(i, P1i, P0i)

        # kann nur am anfang passieren:
        assert has_full_row_rank(P1i), "P10 does not have full row rank. There \
                                        must be algebraic equations."

        # 1. reduktionsschritt
        Bi = reduction(myIteration)

        assert not is_zero_matrix(Bi), "System ist not flat!"

        if outlier(Bi):
            # sonderfall 4.7
           Bi = fourseven(myIteration)

        if end_condition(Bi):
            # flag
            myIteration.last_iter_step = True

            # add to system stack + print
            myStack.add_iteration(myIteration)
            break

        # 2. eliminationsschritt
        P1i, P0i = elimination(myIteration)

        if mode=="manual": print_line()

        # add to system stack + print iteration
        myStack.add_iteration(myIteration)

        i += 1

    # create transformation and calculate Q and G(d/dt)
    myStack.transformation = Transformation(myStack)
    # check for integrability
    myIntegrabilityCheck = IntegrabilityCheck(myStack)

    # for testing
    global P, Q, Q_, G, T, w, I
    T = myStack.transformation
    I = myIntegrabilityCheck
    P = myStack.transformation.P
    Q = myStack.transformation.Q
    # non commutative version:
    Q_ = T.make_symbols_non_commutative(Q)
    G = myStack.transformation.G
    w = T.w

    print_line()

if __name__ == '__main__':
    main()

