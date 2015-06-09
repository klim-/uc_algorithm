# -*- coding: utf-8 -*-

from print_candy import *

from IPython import embed as IPS

class IterationStack(object):
    """ container that carries data for the ith iteration
    """
    def __init__(self, iteration):
        self.i = iteration

        self.is_outlier = False
        self.last_iter_step = False # might not be necessary

        self.P1 = None
        self.P0 = None

        self.P1_roc = None
        self.P1_rpinv = None
        self.P1_dot = None

        self.A = None
        self.B = None
        self.B_loc = None
        self.B_lpinv = None
        self.B_tilde = None

        ## outlier 4.7
        self.Z = None
        self.Z_lpinv = None

    def print_stack(self):
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
            print_matrix("Z", self.i, "", self.Z)
            print_matrix("Z", self.i, "_lpinv", self.Z_lpinv)
