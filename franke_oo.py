# -*- coding: utf-8 -*-

"""
Erste Implementierung des Algorithmus zur Konstruktion flacher Ausgänge für
nichtlineare Systeme von Dr. Matthias Franke, TU Dresden
"""

# enable true divison
from __future__ import division

from IPython import embed as IPS

import numpy as np
import sympy as sp
import symb_tools as st
import diffgeopy as ct


from algebra import *
from examples import *
from integrability import *
from matrix_container import *
from print_candy import *
from system_container import *
from transformations import *


myStack = SystemStack()


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

def reduction(P1i, P0i):

    if roc=="conv":
        P1i_roc = right_ortho_complement(P1i)
    elif roc=="alt":
        P1i_roc = alternative_right_ortho_complement(P1i)

    P1i_rpinv = right_pseudo_inverse(P1i)

    P1i_dot = st.perform_time_derivative(P1i, myStack.diffvec_x)

    Ai = sp.simplify( (P0i - P1i_dot)*P1i_rpinv )
    Bi = sp.simplify( (P0i - P1i_dot)*P1i_roc )

    return Ai, Bi, P1i_roc, P1i_rpinv, P1i_dot

def fourseven(Ai, Bi, P1i_roc, P1i_rpinv):
    if roc=="conv":
        K2 = right_ortho_complement(Bi)
    elif roc=="alt":
        K2 = alternative_right_ortho_complement(Bi)

    if has_full_row_rank(Bi):
        K1 = right_pseudo_inverse(Bi)
    else:
        K1 = right_ortho_complement(K2.T)

    K = st.concat_cols(K1, K2)

    assert is_regular_matrix(K), "K is not a regular matrix."

    Bi_tilde = sp.simplify(Bi*K1) # unit matrix
    P1i_tilde_roc = sp.simplify( P1i_roc*K1 )

    Zi = sp.simplify( P1i_roc*K2 )

    Zi_lpinv = Zi_left_pinv_with_restrictions(Zi, P1i_tilde_roc, P1i_rpinv)

    assert is_unit_matrix( Zi_lpinv*Zi ), "Zi_lpinv seems to be wrong."

    return Zi, Zi_lpinv, Bi_tilde, P1i_tilde_roc

def elimination(Ai, Bi):
    Bi_loc = left_ortho_complement(Bi)
    Bi_lpinv = left_pseudo_inverse(Bi)

    P1i_new = Bi_loc
    P0i_new = Bi_loc*Ai

    return P1i_new, P0i_new, Bi_loc, Bi_lpinv

def store_P_matrices(cls, P1i, P0i):
    cls.P1 = P1i
    cls.P0 = P0i

def store_reduction_matrices(cls, Ai, Bi, P1i_roc, P1i_rpinv, P1i_dot):
    cls.A = Ai
    cls.B = Bi
    cls.P1_roc = P1i_roc
    cls.P1_rpinv = P1i_rpinv
    cls.P1_dot = P1i_dot

def store_outlier_matrices(cls, Zi, Zi_lpinv, Bi_tilde, Pi_tilde_roc):
    cls.Z = Zi
    cls.Z_lpinv = Zi_lpinv
    cls.B_tilde = Bi_tilde
    cls.P1_tilde_roc = Pi_tilde_roc


########################################################################

def main():
    myStack.vec_x = vec_x

    print "x ="; print_nicely(myStack.vec_x)
    print "\n\n xdot ="; print_nicely(myStack.vec_xdot)


    # exterior derivative of F_eq:
    try:
        P1i # in case the system is given by the matrices P1i and P0i
    except NameError:
        print "\n\n0 = F(x,xdot) ="; print_nicely(F_eq)

        P1i = sp.Matrix([])
        P0i = sp.Matrix([])
        for i in xrange(len(myStack.vec_xdot)):
            vector1 = sp.simplify( F_eq.diff(myStack.vec_xdot[i]) )
            P1i = st.concat_cols(P1i, vector1)

            vector0 = sp.simplify( F_eq.diff(myStack.vec_x[i]) )
            P0i = st.concat_cols(P0i, vector0)
    print "\n"*3

    ####################################################################

    i = 0

    while 1:
        # new iteration_stack
        myIteration = IterationStack(i)
        store_P_matrices(myIteration, P1i, P0i)

        # kann nur am anfang passieren:
        assert has_full_row_rank(P1i), "P10 does not have full row rank. There \
                                        must be algebraic equations."

        # 1. reduktionsschritt
        reduction_matrices = reduction(P1i, P0i)
        store_reduction_matrices( myIteration, *reduction_matrices )
        Ai, Bi, P1i_roc, P1i_rpinv, P1i_dot = reduction_matrices

        assert not is_zero_matrix(Bi), "System ist not flat!"

        if outlier(Bi):
            myIteration.is_outlier = True

            outlier_matrices = fourseven(Ai, Bi, P1i_roc, P1i_rpinv)
            store_outlier_matrices( myIteration, *outlier_matrices )
            Zi, Zi_lpinv, Bi_tilde, Pi_tilde_roc = outlier_matrices
            Bi = Bi_tilde

            myIteration.B_tilde_lpinv = left_pseudo_inverse(Bi)
            

        if end_condition(Bi):
            myIteration.last_iter_step = True

            # print results from this iteration
            myIteration.print_stack()

            # add to system stack
            myStack.add_iteration(myIteration)
            break

        # 2. eliminationsschritt
        P1i, P0i, Bi_loc, Bi_lpinv = elimination(Ai, Bi)

        # store
        myIteration.B_loc = Bi_loc
        myIteration.B_lpinv = Bi_lpinv

        # print results from this iteration
        myIteration.print_stack()

        # add to system stack
        myStack.add_iteration(myIteration)

        i += 1

    # create transformation and calculate Q
    # TODO: berechne G(d/dt)
    myStack.transformation = Transformation(myStack)
    # check for integrability
    myIntegrabilityCheck = IntegrabilityCheck(myStack)

    #try_integrability_conditions(w)
    myIntegrabilityCheck.integrability_conditions()


    print_line()

if __name__ == '__main__':
    main()

