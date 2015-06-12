# -*- coding: utf-8 -*-

import sympy as sp
import numpy as np
import symb_tools as st
import non_commutative_tools as nct

from print_candy import *
from algebra import *

from IPython import embed as IPS

# needed noncommutative for G matrix calculation
#t = nct.t
s = sp.Symbol("s", commutative=False)

class Transformation(object):

    def __init__(self, myStack):
        self._myStack = myStack

        self.P = None # P[d/dt]=P0 + P1*(d/dt)

        self.Q = None # unimodular completion of P[d/dt]
        self._w = None # internal variable of corresponding 1form

        self.G = None # G[d/dt]

        self.calculate_P_Matrix()
        self.calculate_Q_matrix()
        self.calculate_G_matrix()

    @property
    def w(self):
        return self._w

    @w.setter
    def w(self, value):
        self._w = value

        # print when set
        print_equation_from_list(self.w, "w", sympy=False)
        print_line()

    def calculate_P_Matrix(self):
        P10 = self._myStack.get_iteration(0).P1
        P00 = self._myStack.get_iteration(0).P0
        
        new_vector = self.get_nc_subs_vector([P10, P00])

        new_vector_nc, subs_tuples = self.make_symbols_non_comutative(new_vector)

        P10_nc = self._myStack.get_iteration(0).P1.subs(self._myStack.diffvec_x, new_vector_nc)
        P00_nc = self._myStack.get_iteration(0).P0.subs(self._myStack.diffvec_x, new_vector_nc)

        self.P = P00_nc + P10_nc*s

    def calculate_Q_matrix(self):
        print_line()
        print("Algorithmus am Ende\n")

        Q_relevant_matrices, Q_tilde_relevant_matrices = self._myStack.get_Q_relevant_matrices()

        def multiply_matrices_in_list(matrix_list):
            result = matrix_list[0]
            for i in xrange(1, len(matrix_list)):
                result = matrix_list[i]*result
            return result

        Q_matrix = multiply_matrices_in_list( Q_relevant_matrices )

        if not Q_tilde_relevant_matrices==[]:
            Q_tilde = multiply_matrices_in_list( Q_tilde_relevant_matrices )
            Q_matrix = st.concat_rows(Q_matrix, Q_tilde)

        m, n = Q_matrix.shape
        assert n==len(self._myStack.vec_x), "Dimensions of Q-Matrix do not fit."

        print "Q-matrix = "
        print_nicely(Q_matrix)
        print "\n"

        self.Q = Q_matrix

    def calculate_G_matrix(self):
        Gi_list = []

        for i in xrange(self._myStack.iteration_steps()):
            Gi = self.calculate_Gi_matrix(i)
            Gi_list.append(Gi)

        G = Gi_list[0]
        for i in xrange(1, len(Gi_list)):
            #G = nct.nc_mul(G,Gi_list[i])
            G = G*Gi_list[i]

        #G = sp.simplify(G)

        # right shift
        

        print "G-matrix = "
        print_nicely(G)

        self.G = G

    def make_symbols_non_comutative(self, expr, appendix='_nc'):
        symbs = st.atoms(expr, sp.Symbol)
        nc_symbols = [s for s in symbs if s.is_commutative]

        new_symbols = [sp.Symbol(s.name+appendix, commutative=False)
                       for s in nc_symbols]

        tup_list = zip(new_symbols, nc_symbols)
        return expr.subs(zip(nc_symbols, new_symbols)), tup_list

    def are_all_symbols_nc(self, matrix):
        for i in list(matrix.atoms(sp.Symbol)):
            if i.is_commutative:
                return False
        return True

    def get_nc_subs_vector(self, matrix_list):
        # get symbols:
        # duplicate from integrability.py (def generate_basis(self))
        # check highest order of derivatives in matrix_list
        highest_order = 0
        for i in xrange(len(matrix_list)):
            for n in matrix_list[i].atoms():
                if hasattr(n, "difforder"):
                    if n.difforder>highest_order:
                        highest_order = n.difforder

        # generate vector with vec_x and its derivatives up to highest_order
        new_vector = self._myStack.vec_x
        for index in xrange(1, highest_order+1):
            vec_x_ndot = st.perform_time_derivative(self._myStack.vec_x, self._myStack.vec_x, order=index)
            new_vector = st.concat_rows( new_vector, vec_x_ndot )

        return new_vector

    def calculate_Gi_matrix(self, i):
        # 1) richtige matrizen für Gi[d/dt] raussuchen
        # 2) symbole durch nichtkommutative ersetzen, Gi[d/dt] ausrechnen
        # 3) G[d/dt] = G0[d/dt] * G1[d/dt] * ... * Gi[d/dt]
        # 4) sonderfallbehandlung 4.7

        iteration = self._myStack.get_iteration(i)
        
        P1_rpinv = iteration.P1_rpinv
        P1_roc = iteration.P1_roc
        B_lpinv = iteration.B_lpinv
        A = iteration.A

        if iteration.is_outlier:
            B_lpinv = iteration.B_tilde_lpinv

        new_vector = self.get_nc_subs_vector([P1_rpinv, P1_roc, B_lpinv, A])

        new_vector_nc, subs_tuples = self.make_symbols_non_comutative(new_vector)

        #IPS()
        # TODO: sonderfallbetrachtung verifizieren

        # substitute symbols
        P1_rpinv_nc = iteration.P1_rpinv.subs(self._myStack.diffvec_x, new_vector_nc)
        A_nc = iteration.A.subs(self._myStack.diffvec_x, new_vector_nc)

        IPS()

        if not iteration.is_outlier:
            P1_roc_nc = iteration.P1_roc.subs(self._myStack.diffvec_x, new_vector_nc)
            B_lpinv_nc = iteration.B_lpinv.subs(self._myStack.diffvec_x, new_vector_nc)

            #Gi = P1_rpinv_nc - P1_roc_nc*( B_lpinv_nc*s + B_lpinv_nc*A_nc )
            Gi = P1_rpinv_nc - nct.nc_mul(P1_roc_nc, nct.nc_mul(B_lpinv_nc,s) + nct.nc_mul(B_lpinv_nc, A_nc))
        else:
            Zi = iteration.Z.subs(self._myStack.diffvec_x, new_vector_nc)

            P1_roc_nc = iteration.P1_tilde_roc.subs(self._myStack.diffvec_x, new_vector_nc)
            B_tilde_lpinv_nc = iteration.B_tilde_lpinv.subs(self._myStack.diffvec_x, new_vector_nc)
            B_lpinv_nc = B_tilde_lpinv_nc

            #Gi = P1_rpinv_nc - P1_roc_nc*( B_lpinv_nc*s + B_lpinv_nc*A_nc )
            Gi = P1_rpinv_nc - nct.nc_mul(P1_roc_nc, nct.nc_mul(B_lpinv_nc,s) + nct.nc_mul(B_lpinv_nc, A_nc))
            Gi = st.concat_cols(Gi, Zi)

        # right shift

        #print_matrix("G", i, "", Gi)

        # store
        iteration.Gi = Gi

        return Gi
















