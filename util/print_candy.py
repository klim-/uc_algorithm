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

def print_special_case_line():
    txt = "--- special case "
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


