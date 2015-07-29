# -*- coding: utf-8 -*-

from print_candy import *

from IPython import embed as IPS

class IterationStack(object):
    """ container that carries data for the ith iteration
    """
    def __init__(self, iteration, P1, P0):
        self.i = iteration

        self.is_outlier = False
        self.last_iter_step = False # might not be necessary

        self.P1 = P1
        self.P0 = P0

        self.P1_roc = None
        self.P1_rpinv = None
        self.P1_dot = None

        self.A = None
        self.B = None
        self.B_loc = None
        self.B_lpinv = None

        ## outlier 4.7
        self.Z = None
        self.Z_lpinv = None
        self.B_tilde = None
        self.B_tilde_lpinv = None
        self.P1_tilde_roc = None

        self.Gi = None

    def store_reduction_matrices(self, A, B, B_lpinv, P1_roc, P1_rpinv, P1_dot):
        self.A = A
        self.B = B
        self.B_lpinv = B_lpinv
        self.P1_roc = P1_roc
        self.P1_rpinv = P1_rpinv
        self.P1_dot = P1_dot

    def store_outlier_matrices(self, Z, Z_lpinv, B_tilde, B_tilde_lpinv, P1_tilde_roc):
        self.is_outlier = True
        self.Z = Z
        self.Z_lpinv = Z_lpinv
        self.B_tilde = B_tilde
        self.B_tilde_lpinv = B_tilde_lpinv
        self.P1_tilde_roc = P1_tilde_roc

    def store_elimination_matrices(self, B_loc):
        self.B_loc = B_loc

    def get_Gi_matrices(self):
        if self.is_outlier:
            return P1_rpinv, P1_tilde_roc, B_tilde_lpinv, A, Z
        else:
            return P1_rpinv, P1_roc, B_lpinv, A

    def print_stack(self):
        """ prints all matrices for this iteration
        """
        print_next_iteration(self.i)

        print_matrix("P1", self.i, "", self.P1)
        print_matrix("P0", self.i, "", self.P0)
        
        print_matrix("P1", self.i, "_roc", self.P1_roc)
        print_matrix("P1", self.i, "_rpinv", self.P1_rpinv)
        print_matrix("P1", self.i, "_dot", self.P1_dot)

        print_matrix("A", self.i, "", self.A)
        print_matrix("B", self.i, "", self.B)
        print_matrix("B", self.i, "_loc", self.B_loc)
        print_matrix("B", self.i, "_lpinv", self.B_lpinv)
        
        if self.is_outlier:
            print_outlier_line()

            print_matrix("B", self.i, "_tilde", self.B_tilde)
            print_matrix("B", self.i, "_tilde_lpinv", self.B_tilde_lpinv)
            print_matrix("P1", self.i, "_tilde_roc", self.P1_tilde_roc)
            print_matrix("Z", self.i, "", self.Z)
            print_matrix("Z", self.i, "_lpinv", self.Z_lpinv)

    def print_hint_stack(self):
        """ prints remaining matrices in case of manual mode
        """
        print_next_iteration(self.i)
        print_matrix("P1", self.i, "_dot", self.P1_dot)

        print_matrix("A", self.i, "", self.A)
        print_matrix("B", self.i, "", self.B)
        print_matrix("B", self.i, "_loc", self.B_loc)
        print_matrix("B", self.i, "_lpinv", self.B_lpinv)
        
        if self.is_outlier:
            print_outlier_line()

            print_matrix("B", self.i, "_tilde", self.B_tilde)
            print_matrix("B", self.i, "_tilde_lpinv", self.B_tilde_lpinv)
            print_matrix("P1", self.i, "_tilde_roc", self.P1_tilde_roc)
            print_matrix("Z", self.i, "", self.Z)
            print_matrix("Z", self.i, "_lpinv", self.Z_lpinv)
