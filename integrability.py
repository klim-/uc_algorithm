# -*- coding: utf-8 -*-

import sympy as sp
import diffgeopy as ct
import symb_tools as st

from print_candy import *

from IPython import embed as IPS

class IntegrabilityCheck(object):
    """ the system is flat iff a one form w exists such that dw=0.
    """
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

    def is_integrable_1(self, w):
        return all(x.d.is_zero() for x in w)

    def is_integrable_2(self, w, i):
        for j in xrange(len(w)):
            if (w[i].d^w[j]).is_zero():
                print "w[" + str(i) + "].d^w[" + str(j) + "]=0"
                if i==j:
                    self.find_pde_w_wedge_wd(w[i])
                return True
            elif (w[i].d^w[j].d).is_zero() and j!=i:
                print "w[" + str(i) + "].d^w[" + str(j) + "].d=0"
                return True
        return False

    def integrability_conditions(self):
        w = self._myStack.transformation.w
        j = len(w)
        y = []
        for i in xrange(j):
            if (w[i].d).is_zero():
                print "w[" + str(i) + "].d = 0"
                print "<-> Integrabilitätsbedingung für w[" + str(i) + "] erfüllt."
                print "    Der Flache Ausgang wurde berechnet:"
                y=0
                for j in xrange(0,len(self._myStack.basis)):
                    y += sp.integrate(w[i].coeff[j], self._myStack.basis[j])
                print "y" + str(i+1) + " ="
                print_nicely(y)
                print "\n\n"
            elif self.is_integrable_2(w, i):
                print "<-> Integrabilitätsbedingung für w[" + str(i) + "] erfüllt.\n\n"
            else:
                print "Integrabilitätsbedingung für w[" + str(i) + "] nicht erfüllt.\n\n"
        
        #if is_integrable_1(w):
        #    integrate_1form(w)

    def integrate_1form(self, w):
        y = []
        for i in xrange(0,len(w)):
            y.append(0)
            for j in xrange(0,len(basis)):
                y[i] += sp.integrate(w[i].coeff[j],basis[j])

        print "\033[93mDer flache Ausgang wurde berechnet:\033[0m\n"
        print_equation_from_list(y, "y")

    def integrate_one_form(self, w_coeff):
        y = 0
        for i in xrange(0,len(w_coeff)):
            y += sp.integrate(w_coeff[i],basis[i])
        
        y=sp.simplify(y)
        
        print "y="
        print_nicely(y)
        return y

    def find_pde_w_wedge_wd(self, w):
        a = w.coeff
        n = len(a) + 1
        pde_list = []
        mu = sp.Function('mu')(*self._myStack.vec_x)
        print "PDGLs wurden berechnet:"
        for i in xrange(1,n-1):
            for j in xrange(i+1,n):
                #print "i=" + str(i) + ", j=" + str(j)
                pde = a[j-1]*mu.diff(self._myStack.vec_x[i-1]) - a[i-1]*mu.diff(self._myStack.vec_x[j-1]) + mu*a[j-1].diff(self._myStack.vec_x[i-1]) - mu*a[i-1].diff(self._myStack.vec_x[j-1])
                #IPS
                print "0 ="
                print_nicely(pde)
                try:
                    sol = sp.pdsolve(pde)
                    print "Integrierender Faktor wurde berechnet:"
                    print_nicely(sol)

                except:
                    print "Could not solve"
                    break
                
                sol = sol.rhs
                
                # berechne einfachsten integrierenden faktor: F(args)=args
                print "Vereinfachter integrierender Faktor:"
                args = sol.args
                funct = None
                for i in xrange(len(args)):
                    if isinstance(args[i], sp.Function):
                        funct = args[i]
                        break
                sol_tilde = sol.subs(funct, *funct.args)
                print_nicely(sol_tilde)
                    

                
                print "Damit ergibt sich die integrierbare Einsform:"
                w_tilde = sol_tilde*w
                print_nicely(w_tilde)
                integrate_one_form(w_tilde.coeff)

                pde_list.append(pde)
        return pde_list

    def reduce_pde(self, pde):
        # works if pde has only one summand
        
        # then the unknown function cannot depend on the differentiation variable
        asd = True
        return True if asd else False

    #def try_integrability_conditions(w):

        ## print output of
        #for i in xrange(len(w)):
            #print "w[" + str(i+1) + "].d = " + str(w[i].d)


        ## dw_i = 0
        #if is_integrable_1(w):
            #print "Integrabilitätsbedingung ist erfüllt <=> dw = 0.\n"
            #y = []
            #for i in xrange(0,len(w)):
                #y.append(0)
                #for j in xrange(0,len(basis)):
                    #y[i] += sp.integrate(w[i].coeff[j],basis[j])

            #print "Der flache Ausgang wurde berechnet:\n"
            #print_equation_from_list(y, "y")
        ## dw_i^w_i = 0
        ##elif:
        ##    pass
            
        #else:
            #print "Integrabilitätsbedingung ist nicht erfüllt <=> dw != 0.\n"

