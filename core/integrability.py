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

import sympy as sp
import diffgeopy as ct
import symb_tools as st

import util.print_candy as pc

from IPython import embed as IPS

class IntegrabilityCheck(object):
    def __init__(self, myStack):
        self._myStack = myStack

        self.generate_basis()
        self.assemble_dual_basis()
        self.integrability_conditions()

    def generate_basis(self):
        # TODO: not the most elegant way?
        # check highest order of derivatives
        highest_order = 0
        for n in self._myStack.transformation.Q.atoms():
            if hasattr(n, "difforder"):
                if n.difforder>highest_order:
                    highest_order = n.difforder

        # generate vector with vec_x and its derivatives up to highest_order
        new_vector = self._myStack.vec_x
        for index in xrange(1, highest_order+1):
            vec_x_ndot = st.perform_time_derivative(self._myStack.vec_x, self._myStack.vec_x, order=index)
            new_vector = st.concat_rows( new_vector, vec_x_ndot )

        # generate basis_1form up to this order
        basis, basis_1form = ct.diffgeo_setup(new_vector)

        # store
        self._myStack.basis = basis
        self._myStack.basis_1form = basis_1form

    def assemble_dual_basis(self):
        _Q = self._myStack.transformation.Q
        p1, p2 = _Q.shape

        l1 = len(self._myStack.basis)

        # w ist a vector, but basis_1form cannot currently be used in sympy
        # matrices, so a list will hold the entries of the column vector w:
        # create list with length p1 (p1 rows)
        w_T = [0] * p1
        for j in xrange(p1):
            # initialize with "empty" 1forms
            w_T[j] = 0 * self._myStack.basis_1form[0]
            for i in xrange(p2):
                w_T[j] = w_T[j] + _Q.row(j)[i] * self._myStack.basis_1form[i]

        self._myStack.transformation.w = w_T

    def is_dwi_zero(self, w, i):
        """ checks if dw_i=0
        """
        return True if (w[i].d).is_zero() else False

    def integrate_dwi(self, w, i):
        print "------"
        print "w[" + str(i) + "].d = 0"
        print "------"
        print "<-> Integrability condition satisfied for w[" + str(i) + "].\n\n"
        print "Integrating w[" + str(i) + "]:"
        y=0
        for j in xrange(0,len(self._myStack.basis)):
            y += sp.integrate(w[i].coeff[j], self._myStack.basis[j])
        print "y[" + str(i) + "] ="
        pc.print_nicely(y)
        print "\n\n"

    def is_dwi_wedge_wi_zero(self, w, i):
        """ checks if dw_i^w_i=0
        """
        return True if (w[i].d^w[i]).is_zero() else False

    def find_pde_dwi_wedge_wi(self, w, i):
        """ calculates pdes dw_i^w_i=0 for some j
        """
        pde_list = [] # a list of pdes that need to be zero

        a = w[i].coeff

        n = len(w[i].basis) + 1
        
        xx = self._myStack.vec_x
        
        mu = sp.Function('mu')(*xx)
        for k in xrange(1, n-1):
            for l in xrange(k+1, n):
                pde =   a[l-1]*mu.diff(xx[k-1]) - a[k-1]*mu.diff(xx[l-1]) +\
                        mu*a[l-1].diff(xx[k-1]) - mu*a[k-1].diff(xx[l-1])
                pde_list.append(pde)
        return pde_list

    def is_dwi_wedge_wi_wedge_wj_zero(self, w, i):
        """ checks if dw_i^w_i^w_j=0 for some j
        """
        for j in xrange(len(w)):
            if (not j==i) and (w[i].d^w[i]^w[j]).is_zero():
                return j
        return False

    def find_pde_dwi_wedge_wi_wedge_wj(self, w, i, j):
        """ calculates pdes dw_i^w_i^w_j=0 for some j
        """
        pde_list = [] # a list of pdes that need to be zero
        
        a = w[i].coeff
        b = w[j].coeff
        
        n = len(w[i].basis) + 1
        
        xx = self._myStack.vec_x
        
        mu1 = sp.Function('mu_1')(*xx)
        mu2 = sp.Function('mu_2')(*xx)
        for k in xrange(1, n-1):
            for l in xrange(k+1, n):
                pde =   a[l-1]*mu1.diff(xx[k-1]) - a[k-1]*mu1.diff(xx[l-1]) +\
                        mu1*a[l-1].diff(xx[k-1]) - mu1*a[k-1].diff(xx[l-1]) +\
                        b[l-1]*mu2.diff(xx[k-1]) - b[k-1]*mu2.diff(xx[l-1]) +\
                        mu2*b[l-1].diff(xx[k-1]) - mu2*b[k-1].diff(xx[l-1])
                pde_list.append(pde)
        return pde_list

    def integrability_conditions(self):
        """ this function checks the following integrability conditions:
                        dw_i=0
                    dw_i^w_i=0
            and dw_i^w_i^w_j=0
        """
        w = self._myStack.transformation.w
        #~ j = len(w)
        for i in xrange(len(w)):
            if self.is_dwi_zero(w, i):
                self.integrate_dwi(w, i)
            elif self.is_dwi_wedge_wi_zero(w,i):
                print "------"
                print "w[" + str(i) + "].d^w[" + str(i) + "] = 0"
                print "------"
                print "resp."
                print "du=0 for u=mu1*w["+str(i)+"]"
                print "<-> Integrability condition satisfied for w[" + str(i) + "].\n\n"
                var = raw_input("Type \"pde\" to calculate PDE-conditions for mu1\n")
                if var=="pde":
                    pde_list = self.find_pde_dwi_wedge_wi(w,i)
                    pc.print_pde_list(pde_list)
            # TODO: this is not very elegant (is_dwi_wedge_wi_wedge_wj_zero
            #       returns integer or False)
            elif self.is_dwi_wedge_wi_wedge_wj_zero(w, i)!=False:
                j = self.is_dwi_wedge_wi_wedge_wj_zero(w, i)
                print "------"
                print "w[" + str(i) + "].d^w[" + str(i) + "]^w[" + str(j) +"] = 0"
                print "------"
                print "resp."
                print "du=0 for u=mu1*w["+str(i)+"] + mu2*w["+str(j)+"]"
                print "<-> Integrability condition satisfied for w[" + str(i) + "].\n\n"
                var = raw_input("Type \"pde\" to calculate PDE-conditions for mu1 and mu2\n")
                if var=="pde":
                    pde_list = self.find_pde_dwi_wedge_wi_wedge_wj(w,i,j)
                    pc.print_pde_list(pde_list)
            else:
                print "------"
                print "Integrability condition not satisfied for w[" + str(i) + "]. Please try manually."
                print "------\n\n"
