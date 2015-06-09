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
from system import *
from system_container import *
from transformations import *


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

def Zi_left_pinv_with_restrictions(Zi, P1i_tilde_right_ortho, P1i_right_pseudo_inv):
    """ Given a matrix Zi, this function calculates a matrix such that:
            Zi_left_pinv * Zi                    = I
            Zi_left_pinv * P1i_tilde_right_ortho = 0
            Zi_left_pinv * P1i_right_pseudo_inv  = 0
    """
    assert check_row_compatibility(Zi, P1i_tilde_right_ortho, P1i_right_pseudo_inv),\
        "Matrices do not match in row dimension."

    C = st.concat_cols(Zi, P1i_tilde_right_ortho, P1i_right_pseudo_inv)

    assert is_regular_matrix(C), "C is not a regular matrix"
    C_det = C.berkowitz_det()
    C_inv = C.adjugate()/C_det
    C_inv = sp.simplify(C.inv())

    m, n = Zi.shape
    Zi_left_pinv = sp.Matrix([])
    for i in xrange(0,n):
        Zi_left_pinv = st.concat_rows(Zi_left_pinv,C_inv.row(i))

    o, p = Zi_left_pinv.shape
    assert o==n and p==m, "There must have been a problem with the\
                            computation of Zi_left_pinv"

    assert is_unit_matrix(Zi_left_pinv*Zi), "Zi_left_pinv is wrong"
    assert is_zero_matrix(Zi_left_pinv*P1i_tilde_right_ortho), "Zi_left_pinv is wrong"
    assert is_zero_matrix(Zi_left_pinv*P1i_right_pseudo_inv), "Zi_left_pinv is wrong"

    return Zi_left_pinv

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

    #K1 = right_ortho_complement(K2.T)
    K1 = right_pseudo_inverse(Bi)

    K = st.concat_cols(K1, K2)

    assert is_regular_matrix(K), "K is not a regular matrix."

    Bi_tilde = sp.simplify(Bi*K1) # unit matrix
    P1i_tilde_roc = sp.simplify( P1i_roc*K1 )
    Zi = sp.simplify( P1i_roc*K2 )

    Zi_lpinv = Zi_left_pinv_with_restrictions(Zi, P1i_tilde_roc, P1i_rpinv)

    assert is_unit_matrix( Zi_lpinv*Zi ), "Zi_lpinv seems to be wrong."

    return Zi, Zi_lpinv, Bi_tilde

def elimination(Ai, Bi):
    Bi_loc = left_ortho_complement(Bi)
    Bi_lpinv = left_pseudo_inverse(Bi)

    P1i_new = Bi_loc
    P0i_new = Bi_loc*Ai

    return P1i_new, P0i_new, Bi_loc, Bi_lpinv

def set_P_matrices(cls, P1i, P0i):
    cls.P1 = P1i
    cls.P0 = P0i

def set_reduction_matrices(cls, Ai, Bi, P1i_roc, P1i_rpinv, P1i_dot):
    cls.A = Ai
    cls.B = Bi
    cls.P1_roc = P1i_roc
    cls.P1_rpinv = P1i_rpinv
    cls.P1_dot = P1i_dot

def set_outlier_matrices(cls, Zi, Zi_lpinv, Bi_tilde):
    cls.Z = Zi
    cls.Z_lpinv = Zi_lpinv
    cls.B_tilde = Bi_tilde


########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################

#vec_xdot = st.perform_time_derivative(vec_x,vec_x)
#diffvec_x = st.concat_rows(vec_x, vec_xdot)

def main():
    global myStack
    myStack = SystemStack()
    
    myStack.vec_x = vec_x

    print "x ="; print_nicely(myStack.vec_x)
    print "\n\n xdot ="; print_nicely(myStack.vec_xdot)


    # exterior derivative of F_eq:
    # TODO: Not elegant, could be combined!
    try:
        P1i # in case the system is given by the matrices P1i and P0i
    except NameError:
        print "\n\n0 = F(x,xdot) ="; print_nicely(F_eq)
        P1i = sp.Matrix([])
        for i in xrange(len(myStack.vec_xdot)):
            vector = sp.simplify( F_eq.diff(myStack.vec_xdot[i]) )
            P1i = st.concat_cols(P1i, vector)

        P0i = sp.Matrix([])
        for i in xrange(len(myStack.vec_x)):
            vector = sp.simplify( F_eq.diff(myStack.vec_x[i]) )
            P0i = st.concat_cols(P0i, vector)
    print "\n\n\n"

    ####################################################################
    ####################################################################
    ####################################################################
    ####################################################################

    i = 0


    while 1:
        # new iteration_stack
        myIteration = IterationStack(i)
        set_P_matrices(myIteration, P1i, P0i)

        #print_matrix("while P0",i,"",P0i)
        # kann nur am anfang passieren:
        assert has_full_row_rank(P1i), "P10 does not have full row rank. There \
                                        must be algebraic equations."

        # 1. reduktionsschritt
        reduction_matrices = reduction(P1i, P0i)
        set_reduction_matrices( myIteration, *reduction_matrices )
        Ai, Bi, P1i_roc, P1i_rpinv, P1i_dot = reduction_matrices
        

        assert not is_zero_matrix(Bi), "System ist not flat!"

        if outlier(Bi):
            myIteration.is_outlier = True

            outlier_matrices = fourseven(Ai, Bi, P1i_roc, P1i_rpinv)
            set_outlier_matrices( myIteration, *outlier_matrices )
            Zi, Zi_lpinv, Bi_tilde = outlier_matrices
            Bi = Bi_tilde

        if end_condition(Bi):
            myIteration.print_stack()

            myIteration.last_iter_step = True
            myStack.add_iteration(myIteration)
            break

        # 2. eliminationsschritt
        P1i, P0i, Bi_loc, Bi_lpinv = elimination(Ai, Bi)
        myIteration.B_loc = Bi_loc
        myIteration.B_lpinv = Bi_lpinv

        myIteration.print_stack()
        myStack.add_iteration(myIteration)

        i += 1

    # berechne neues Q
    Q = calculate_Q_matrix(myStack)
    
    #IPS()

    # TODO: not the most elegant way?
    # check highest order of derivatives
    highest_order = 0
    for n in Q.atoms():
        if hasattr(n, "difforder"):
            if n.difforder>highest_order:
                highest_order = n.difforder

    # generate vector with vec_x and its derivatives up to highest_order
    new_vector = myStack.vec_x
    for index in xrange(1, highest_order+1):
        vec_x_ndot = st.perform_time_derivative(myStack.vec_x, myStack.vec_x, order=index)
        new_vector = st.concat_rows( new_vector, vec_x_ndot )

    # generate basis_1form up to this order
    # TODO: there must be a more elegant way for this!
    global basis
    global basis_1form
    basis, basis_1form = ct.diffgeo_setup(new_vector)

    # calculate w
    global w
    w = assemble_dual_basis(Q, basis, basis_1form)

    print_equation_from_list(w, "w", sympy=False)
    print_line()

    #try_integrability_conditions(w)
    integrability_conditions(w, basis)

    print_line()

if __name__ == '__main__':
    main()

