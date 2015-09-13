# -*- coding: utf-8 -*-

#    Copyright (C) 2015
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

import util.print_candy as pc

from IPython import embed as IPS

class IterationStack(object):
    """ container that carries data for the ith iteration
    """
    def __init__(self, iteration, P1, P0):
        self.i = iteration

        self.is_special_case = False
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

        # special case
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

    def store_special_case_matrices(self, Z, Z_lpinv, B_tilde, B_tilde_lpinv, P1_tilde_roc):
        self.is_special_case = True
        self.Z = Z
        self.Z_lpinv = Z_lpinv
        self.B_tilde = B_tilde
        self.B_tilde_lpinv = B_tilde_lpinv
        self.P1_tilde_roc = P1_tilde_roc

    def store_elimination_matrices(self, B_loc):
        self.B_loc = B_loc

    def get_Gi_matrices(self):
        if self.is_special_case:
            return P1_rpinv, P1_tilde_roc, B_tilde_lpinv, A, Z
        else:
            return P1_rpinv, P1_roc, B_lpinv, A

    def print_stack(self):
        """ prints all matrices for this iteration
        """
        pc.print_next_iteration(self.i)

        pc.print_matrix("P1", self.i, "", self.P1)
        pc.print_matrix("P0", self.i, "", self.P0)
        
        pc.print_matrix("P1", self.i, "_roc", self.P1_roc)
        pc.print_matrix("P1", self.i, "_rpinv", self.P1_rpinv)
        pc.print_matrix("P1", self.i, "_dot", self.P1_dot)

        pc.print_matrix("A", self.i, "", self.A)
        pc.print_matrix("B", self.i, "", self.B)
        pc.print_matrix("B", self.i, "_loc", self.B_loc)
        pc.print_matrix("B", self.i, "_lpinv", self.B_lpinv)
        
        if self.is_special_case:
            pc.print_special_case_line()

            pc.print_matrix("B", self.i, "_tilde", self.B_tilde)
            pc.print_matrix("B", self.i, "_tilde_lpinv", self.B_tilde_lpinv)
            pc.print_matrix("P1", self.i, "_tilde_roc", self.P1_tilde_roc)
            pc.print_matrix("Z", self.i, "", self.Z)
            pc.print_matrix("Z", self.i, "_lpinv", self.Z_lpinv)

    def print_hint_stack(self):
        """ prints remaining matrices in case of manual mode
        """
        pc.print_next_iteration(self.i)
        pc.print_matrix("P1", self.i, "_dot", self.P1_dot)

        pc.print_matrix("A", self.i, "", self.A)
        pc.print_matrix("B", self.i, "", self.B)
        pc.print_matrix("B", self.i, "_loc", self.B_loc)
        pc.print_matrix("B", self.i, "_lpinv", self.B_lpinv)
        
        if self.is_special_case:
            pc.print_special_case_line()

            pc.print_matrix("B", self.i, "_tilde", self.B_tilde)
            pc.print_matrix("B", self.i, "_tilde_lpinv", self.B_tilde_lpinv)
            pc.print_matrix("P1", self.i, "_tilde_roc", self.P1_tilde_roc)
            pc.print_matrix("Z", self.i, "", self.Z)
            pc.print_matrix("Z", self.i, "_lpinv", self.Z_lpinv)
