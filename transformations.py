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
        var = raw_input("Type \"G\" to calculate G(d/dt)-matrix\n")
        if var=="G": self.calculate_G_matrix()

    @property
    def w(self):
        return self._w

    @w.setter
    def w(self, value):
        self._w = value

        # print when set
        print_equation_from_list(self.w, "w", sympy=False)
        print_line()

    def mul_rs_and_simplify_nc(self, A, B):
        """ multiplies nc matrices A and B, right shifts s, converts
            the result to commutative and simplifies it
        """
        AB_nc = A*B
        AB_nc_shifted = self.right_shift_all_in_matrix(AB_nc)
        AB_shifted_c = self.make_symbols_commutative(AB_nc_shifted)
        extr_result = list(AB_shifted_c)[0]
        result = custom_simplify(extr_result)

        return result

    def make_symbols_non_commutative(self, expr):
        """ replaces commutative symbols x with non commutative symbols z
        """
        symbs = st.atoms(expr, sp.Symbol)
        nc_symbols = [s for s in symbs if s.is_commutative]

        new_symbols = [sp.Symbol(s.name.replace("x","z"), commutative=False)
                       for s in nc_symbols]

        tup_list = zip(nc_symbols, new_symbols)
        return expr.subs(zip(nc_symbols, new_symbols))#, tup_list

    def convert_to_nc_matrices(self, *args):
        """ Converts commutative symbols x with non commutative symbols z
            in matrices.
        """
        # TODO: asserts!
        
        list_of_nc_matrices = []
        
        for a in args:
            m = self.make_symbols_non_commutative(a)
            list_of_nc_matrices.append(m)
        
        return tuple(list_of_nc_matrices)

    def make_symbols_commutative(self, expr):
        """
        :param expr:
        :return: expr (with all symbols commutative) and
                  a subs_tuple_list [(s1_c, s1_nc), ... ]
        """

        symbs = st.atoms(expr, sp.Symbol)
        nc_symbols = [s for s in symbs if not s.is_commutative]

        new_symbols = [sp.Symbol(s.name.replace("z","x"), commutative=True)
                       for s in nc_symbols]

        tup_list = zip(new_symbols, nc_symbols)
        return expr.subs(zip(nc_symbols, new_symbols))#, tup_list

    def all_nc(self, matrix):
        """ checks if matrix is fully non commutative
        """
        for i in list(matrix.atoms(sp.Symbol)):
            if i.is_commutative:
                return False
        # wenn es eine reine Zahlenmatrix ist kommt ebenfalls True raus
        # ist das gut?
        return True

    def right_shift_all_in_matrix(self, matrix):
        # @carsten: ist diese funktion wirklich notwendig?
        m,n = matrix.shape
        matrix_shifted = sp.Matrix([])
        t_dep_symbols = [symb for symb in st.atoms(matrix, sp.Symbol) if not symb == s]
        for i in xrange(n):
            col = matrix.col(i)
            col_shifted = sp.Matrix([nct.right_shift_all(expr,s,t, t_dep_symbols) for expr in col])
            matrix_shifted = st.concat_cols(matrix_shifted, col_shifted)

        return matrix_shifted

########################################################################

    def calculate_P_Matrix(self):
        """ specifies system with P(d/dt) = P00 + P10*s
        """
        P10 = self._myStack.get_iteration(0).P1
        P00 = self._myStack.get_iteration(0).P0

        P10_nc = self.make_symbols_non_commutative(P10)
        P00_nc = self.make_symbols_non_commutative(P00)

        self.P = P00_nc + nct.nc_mul(P10_nc, s)
        return self.P

    def multiply_matrices_in_list(self, matrix_list):
        result = matrix_list[0]
        for i in xrange(1, len(matrix_list)):
            result = matrix_list[i]*result
        return result

    def calculate_Q_matrix(self):
        """ unimodular completion of P(d/dt) with Q = P1i * ... * P11 * P10
        """
        print_line()
        print("Algorithmus am Ende\n")

        Q_relevant_matrices, Q_tilde_relevant_matrices = self._myStack.get_Q_relevant_matrices()


        #Q_matrix = self.multiply_matrices_in_list( Q_relevant_matrices )

        #if not len(Q_tilde_relevant_matrices)==0:
            #Q_tilde = self.multiply_matrices_in_list( Q_tilde_relevant_matrices )
            #Q_matrix = st.concat_rows(Q_matrix, Q_tilde)
        Q1 = self.multiply_matrices_in_list( Q_relevant_matrices )

        if not len(Q_tilde_relevant_matrices)==0:
            Q2 = self.multiply_matrices_in_list( Q_tilde_relevant_matrices )
        else:
            Q2 = sp.Matrix([])

        Q = st.concat_rows(Q1, Q2)

        m, n = Q.shape
        assert n==len(self._myStack.vec_x), "Dimensions of Q-Matrix do not fit."

        print "Q-matrix = "
        print_nicely(Q)
        print "\n"

        self.Q = Q

    def calculate_G_matrix(self):
        G = self.calculate_Gi_matrix(0)
        for i in xrange(self._myStack.iteration_steps()-1):
            G = G*self.calculate_Gi_matrix(i+1)

        # right shift
        G_shifted = self.right_shift_all_in_matrix(G)

        print "G-matrix = "
        print_nicely(G_shifted)

        self.G = G_shifted

    def calculate_Gi_matrix(self, i):
        # 1) richtige matrizen für Gi[d/dt] raussuchen
        # 2) symbole durch nichtkommutative ersetzen, Gi[d/dt] ausrechnen
        # 3) G[d/dt] = G0[d/dt] * G1[d/dt] * ... * Gi[d/dt]
        # 4) sonderfallbehandlung 4.7 (ähnlich)
        iteration = self._myStack.get_iteration(i)

        if not iteration.is_outlier:
            # get matrices
            P1_rpinv = iteration.P1_rpinv
            P1_roc = iteration.P1_roc
            B_lpinv = iteration.B_lpinv
            A = iteration.A

            # convert commutative symbols to non commutative
            #P1_rpinv_nc = list( self.make_symbols_non_commutative(P1_rpinv) )[0]
            
            P1_rpinv_nc, P1_roc_nc, B_lpinv_nc, A_nc = self.convert_to_nc_matrices(P1_rpinv, P1_roc, B_lpinv, A)
            
            #~ P1_rpinv_nc = self.make_symbols_non_commutative(P1_rpinv)
            #~ P1_roc_nc = self.make_symbols_non_commutative(P1_roc)
            #~ B_lpinv_nc = self.make_symbols_non_commutative(B_lpinv)
            #~ A_nc = self.make_symbols_non_commutative(A)

            Gi = P1_rpinv_nc - P1_roc_nc*( B_lpinv_nc*s_ + B_lpinv_nc*A_nc )
            #Gi = P1_rpinv_nc - nct.nc_mul(P1_roc_nc, nct.nc_mul(B_lpinv_nc,s) + nct.nc_mul(B_lpinv_nc, A_nc))

        else:
            # get matrices
            P1_rpinv = iteration.P1_rpinv
            P1_tilde_roc = iteration.P1_tilde_roc
            B_tilde_lpinv = iteration.B_tilde_lpinv
            A = iteration.A
            Z = iteration.Z

            # convert commutative symbols to non commutative
            P1_rpinv_nc = self.make_symbols_non_commutative(P1_rpinv)
            P1_tilde_roc_nc = self.make_symbols_non_commutative(P1_tilde_roc)
            B_tilde_lpinv_nc = self.make_symbols_non_commutative(B_tilde_lpinv)
            A_nc = self.make_symbols_non_commutative(A)
            Z_nc = self.make_symbols_non_commutative(Z)

            Gi = P1_rpinv_nc - P1_tilde_roc_nc*( B_tilde_lpinv_nc*s_ + B_tilde_lpinv_nc*A_nc )
            #Gi = P1_rpinv_nc - nct.nc_mul(P1_tilde_roc_nc, nct.nc_mul(B_tilde_lpinv_nc,s) + nct.nc_mul(B_tilde_lpinv_nc, A_nc))

            Gi = st.concat_cols(Gi, Z_nc)

        if False: #show_Gi_matrices:
            print_matrix("G", i, "", Gi)

        # store
        iteration.Gi = Gi

        return Gi
