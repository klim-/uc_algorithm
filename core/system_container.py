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

import numpy as np
import sympy as sp
import symb_tools as st
import pycartan as ct

import core.matrix_container as mc

from IPython import embed as IPS

class SystemStack(object):
    """ container that carries resulting matrizes for the system.
        no calculations here, just store and access of data
    """
    def __init__(self, diff_symbols):
        self.iteration_data = [] # contains a matrix container for every iteration step (=index)

        self._vec_x = None # internal variable!

        self.vec_xdot = None
        self.diffvec_x = None # consists of vec_x and vec_xdot
        self.diff_symbols = diff_symbols

        # basis and basis 1-form
        self.basis = None
        self.basis_1form = None

        # transformations
        self.transformation = None

        # quick hack (should probably not live here)
        subs_back_to_c = []

    @property
    def vec_x(self):
        return self._vec_x

    @vec_x.setter
    def vec_x(self, value):
        self._vec_x = value

        # calculate vec_xdot
        self.vec_xdot = st.perform_time_derivative(self.vec_x, self.vec_x)

        # vector for differentiation
        self.diffvec_x = st.concat_rows(self.vec_x, self.vec_xdot, self.diff_symbols)

    def add_iteration(self, iter_stack):
        assert isinstance(iter_stack, mc.IterationStack)
        self.iteration_data.append(iter_stack)
        iter_stack.print_stack()

    def get_iteration(self, i):
        assert type(i)==int, "Iteration step i must be of type integer!"
        return self.iteration_data[i]

    def iteration_steps(self):
        return len(self.iteration_data)

    def print_iteration(self, i):
        self.iteration_data[i].print_stack()

    def get_special_cases(self):
        #~ counter = 0
        special_cases = [0]*self.iteration_steps() # list [0,0,0,...,1,..,0] (1:=special_case)
        for i in xrange(self.iteration_steps()):
            if self.iteration_data[i].is_special_case == True:
                #~ counter += 1
                special_cases[i] = 1
        return special_cases

    def special_case_in_between(self, special_cases):
        if (special_cases.count(1)>0) and (special_cases.index(1)<(len(special_cases))):
            return True
        else:
            return False

    def get_Q_relevant_matrices(self):
        """ all P1i and Zi_lpinv matrices in correct order, so these just have to
            be multiplied to get Q
        """
        Q_relevant_matrices = []
        for i in xrange(self.iteration_steps()):
            P1i = self.iteration_data[i].P1
            Q_relevant_matrices.append(P1i)

        Q_tilde_relevant_matrices = []

        special_cases = self.get_special_cases()

        if special_cases.count(1) == 1:
            special_cases_iteration = special_cases.index(1)
            for i in xrange(self.iteration_steps()):
                if i==special_cases_iteration:
                    Zi_lpinv = self.iteration_data[i].Z_lpinv
                    Q_tilde_relevant_matrices.append(Zi_lpinv)
                else:
                    P1i = self.iteration_data[i].P1
                    Q_tilde_relevant_matrices.append(P1i)
        elif self.special_case_in_between(special_cases):
            raise NotImplementedError

        return Q_relevant_matrices, Q_tilde_relevant_matrices
