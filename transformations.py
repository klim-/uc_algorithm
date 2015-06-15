# -*- coding: utf-8 -*-

import sympy as sp
import numpy as np
import symb_tools as st
import non_commutative_tools as nct

from print_candy import *
from algebra import *

from IPython import embed as IPS

# needed noncommutative for G matrix calculation
t = nct.t
#s_ = nct.s
s = sp.Symbol("s", commutative=False)
s_=s
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

        # test:
        self.PG_nc_shifted = self.mul_and_rs(self.P, self.G)

        # substitute back
        #self.PG = nct.make_all_symbols_commutative(self.PG_nc_shifted)[0]

    @property
    def w(self):
        return self._w

    @w.setter
    def w(self, value):
        self._w = value

        # print when set
        print_equation_from_list(self.w, "w", sympy=False)
        print_line()

    def mul_and_rs(self, A, B):
        # just for quicker testing
        AB = nct.nc_mul( A, B )
        return self.right_shift_all_in_matrix(AB)

    def calculate_P_Matrix(self):
        P10 = self._myStack.get_iteration(0).P1
        P00 = self._myStack.get_iteration(0).P0
        
        new_vector = self.get_nc_subs_vector([P10, P00])
        new_vector_nc, subs_tuples = self.make_symbols_non_comutative(new_vector)

        P10_nc = st.symbs_to_func( P10.subs(subs_tuples), new_vector_nc, t )
        P00_nc = st.symbs_to_func( P00.subs(subs_tuples), new_vector_nc, t )

        #P10_nc = self._myStack.get_iteration(0).P1.subs(self._myStack.diffvec_x, new_vector_nc)
        #P00_nc = self._myStack.get_iteration(0).P0.subs(self._myStack.diffvec_x, new_vector_nc)

        self.P = P00_nc + nct.nc_mul(P10_nc, s_)

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
            G = G*Gi_list[i]
            #G = nct.nc_mul(G, Gi_list[i])

        #G = sp.simplify(G)

        # right shift
        G_shifted = self.right_shift_all_in_matrix(G)
        

        print "G-matrix = "
        print_nicely(G_shifted)

        self.G = G_shifted

    def make_symbols_non_comutative(self, expr, appendix='_nc'):
        symbs = st.atoms(expr, sp.Symbol)
        nc_symbols = [s for s in symbs if s.is_commutative]

        new_symbols = [sp.Symbol(s.name+appendix, commutative=False)
                       for s in nc_symbols]

        tup_list = zip(nc_symbols, new_symbols)
        return expr.subs(zip(nc_symbols, new_symbols)), tup_list

    def all_nc(self, matrix):
        for i in list(matrix.atoms(sp.Symbol)):
            if i.is_commutative:
                return False
        # wenn es eine reine Zahlenmatrix ist kommt ebenfalls True raus
        # ist das gut?
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

        # hack:
        highest_order = 5

        # generate vector with vec_x and its derivatives up to highest_order
        new_vector = self._myStack.vec_x
        for index in xrange(1, highest_order+1):
            vec_x_ndot = st.perform_time_derivative(self._myStack.vec_x, self._myStack.vec_x, order=index)
            new_vector = st.concat_rows( new_vector, vec_x_ndot )

        return new_vector

    def right_shift_all_in_matrix(self, matrix):
        # @carsten: ist diese funktion wirklich notwendig?
        m,n = matrix.shape
        # TODO: im moment nur spaltenmatrix
        matrix_shifted = sp.Matrix([])
        for i in xrange(m):
            row = matrix.row(i)
            row_shifted = sp.Matrix([])
            for j in xrange(n):
                #IPS()
                #sp.pprint(row[j])
                new_row = nct.right_shift_all(row[j], s_, t)
                row_shifted = st.concat_cols(row_shifted, new_row)
            matrix_shifted = st.concat_rows(matrix_shifted, row_shifted)


        #for i in xrange(len(matrix)):
            #new_row = nct.right_shift_all( matrix[i], s_, t )
            #matrix_shifted = st.concat_rows( matrix_shifted, new_row )
        #IPS()
        return matrix_shifted

    def calculate_Gi_matrix(self, i):
        # 1) richtige matrizen für Gi[d/dt] raussuchen
        # 2) symbole durch nichtkommutative ersetzen, Gi[d/dt] ausrechnen
        # 3) G[d/dt] = G0[d/dt] * G1[d/dt] * ... * Gi[d/dt]
        # 4) sonderfallbehandlung 4.7 (ähnlich)
        iteration = self._myStack.get_iteration(i)

        if not iteration.is_outlier:
            P1_rpinv = iteration.P1_rpinv
            P1_roc = iteration.P1_roc
            B_lpinv = iteration.B_lpinv
            A = iteration.A

            new_vector = self.get_nc_subs_vector([P1_rpinv, P1_roc, B_lpinv, A])
            new_vector_nc, subs_tuples = self.make_symbols_non_comutative(new_vector)
           
            P1_rpinv_nc = st.symbs_to_func( P1_rpinv.subs(subs_tuples), new_vector_nc, t )
            P1_roc_nc = st.symbs_to_func( P1_roc.subs(subs_tuples), new_vector_nc, t )
            B_lpinv_nc = st.symbs_to_func( B_lpinv.subs(subs_tuples), new_vector_nc, t )
            A_nc = st.symbs_to_func( A.subs(subs_tuples), new_vector_nc, t )

            #IPS()

            Gi = P1_rpinv_nc - P1_roc_nc*( B_lpinv_nc*s_ + B_lpinv_nc*A_nc )
            #Gi = P1_rpinv_nc - nct.nc_mul(P1_roc_nc, nct.nc_mul(B_lpinv_nc,s) + nct.nc_mul(B_lpinv_nc, A_nc))

        else:
            P1_rpinv = iteration.P1_rpinv
            P1_tilde_roc = iteration.P1_tilde_roc
            B_tilde_lpinv = iteration.B_tilde_lpinv
            A = iteration.A

            Z = iteration.Z

            new_vector = self.get_nc_subs_vector([P1_rpinv, P1_tilde_roc, B_tilde_lpinv, A, Z])
            new_vector_nc, subs_tuples = self.make_symbols_non_comutative(new_vector)
            
            P1_rpinv_nc = st.symbs_to_func( P1_rpinv.subs(subs_tuples), new_vector_nc, t )
            P1_tilde_roc_nc = st.symbs_to_func( P1_tilde_roc.subs(subs_tuples), new_vector_nc, t )
            B_tilde_lpinv_nc = st.symbs_to_func( B_tilde_lpinv.subs(subs_tuples), new_vector_nc, t )
            A_nc = st.symbs_to_func( A.subs(subs_tuples), new_vector_nc, t )
            
            Z_nc = st.symbs_to_func( Z.subs(subs_tuples), new_vector_nc, t )

            #IPS()


            Gi = P1_rpinv_nc - P1_tilde_roc_nc*( B_tilde_lpinv_nc*s_ + B_tilde_lpinv_nc*A_nc )
            #Gi = P1_rpinv_nc - nct.nc_mul(P1_tilde_roc_nc, nct.nc_mul(B_tilde_lpinv_nc,s) + nct.nc_mul(B_tilde_lpinv_nc, A_nc))

            Gi = st.concat_cols(Gi, Z_nc)

        print_matrix("G", i, "", Gi)

        # store
        iteration.Gi = Gi

        return Gi

################
        
        #P1_rpinv = iteration.P1_rpinv
        #P1_roc = iteration.P1_roc
        #B_lpinv = iteration.B_lpinv
        #A = iteration.A

        #if iteration.is_outlier:
            #B_lpinv = iteration.B_tilde_lpinv

        #new_vector = self.get_nc_subs_vector([P1_rpinv, P1_roc, B_lpinv, A])
        #new_vector_nc, subs_tuples = self.make_symbols_non_comutative(new_vector)

        ##IPS()

        ## substitute symbols
        #P1_rpinv_nc = P1_rpinv.subs(subs_tuples)
        #A_nc = A.subs(subs_tuples)

        #IPS()

        #if not iteration.is_outlier:
            #P1_roc_nc = iteration.P1_roc.subs(subs_tuples)
            #B_lpinv_nc = iteration.B_lpinv.subs(subs_tuples)

            #Gi = P1_rpinv_nc - P1_roc_nc*( B_lpinv_nc*s + B_lpinv_nc*A_nc )
            ##Gi = P1_rpinv_nc - nct.nc_mul(P1_roc_nc, nct.nc_mul(B_lpinv_nc,s) + nct.nc_mul(B_lpinv_nc, A_nc))
        #else:
            #Zi = iteration.Z.subs(subs_tuples)

            #P1_roc_nc = iteration.P1_tilde_roc.subs(subs_tuples)
            #B_tilde_lpinv_nc = iteration.B_tilde_lpinv.subs(subs_tuples)
            #B_lpinv_nc = B_tilde_lpinv_nc

            #Gi = P1_rpinv_nc - P1_roc_nc*( B_lpinv_nc*s + B_lpinv_nc*A_nc )
            ##Gi = P1_rpinv_nc - nct.nc_mul(P1_roc_nc, nct.nc_mul(B_lpinv_nc,s) + nct.nc_mul(B_lpinv_nc, A_nc))
            #Gi = st.concat_cols(Gi, Zi)

        ## right shift

        #print_matrix("G", i, "", Gi)

        ## store
        #iteration.Gi = Gi

        #return Gi

















