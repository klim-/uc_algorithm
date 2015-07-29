# -*- coding: utf-8 -*-

import sympy as sp

# enable latex printing in ipython qtconsole
qt = 1
if qt:
    from sympy import init_printing
    from IPython.display import display
    init_printing()

string_line = "#"*85 + "\n\n"
string_line2 = "-"*85 + "\n\n"

def print_matrix(name, i, supplement, matrix):
    if not matrix==None:
        print name + str(i) + supplement + " [" + str(matrix.shape[0]) + " x "  + str(matrix.shape[1]) + "] = "
        print_nicely(matrix)

def print_nicely(formula):
    if qt:
        display(formula)
    else:
        sp.pprint(formula)
    print("\n")

def print_next_iteration(i):
    print "i = " + str(i) + " " + string_line[len(str(i))+5:]

def print_outlier_line():
    txt = "--- Sonderfall 4.7 "
    print txt + string_line2[len(txt):]

def print_line():
    print string_line

def print_equation_from_list(f, var, sympy=True):
    if sympy==False:
        for k in xrange(0, len(f)):
            print str(var) + "[" + str(k) + "] = " + str(f[k]) + "\n"
    else:
        for k in xrange(0, len(f)):
            print str(var) + "[" + str(k) + "] ="
            print_nicely(f[k])

def print_pde_list(pde_list, sympy=True):
    if sympy==False:
        for k in xrange(0, len(pde_list)):
            print "0 = " + str(pde_list[k]) + "\n"
    else:
        for k in xrange(0, len(pde_list)):
            print "0 = "
            print_nicely(pde_list[k])


