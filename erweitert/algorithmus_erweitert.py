# -*- coding: utf-8 -*-

"""
Erste Implementierung des erweiterten Algorithmus zur Konstruktion flacher Ausgänge für
nichtlineare Systeme von Dr. Matthias Franke, TU Dresden
"""

# enable true divison
from __future__ import division

import sympy as sp
import symb_tools as st
import diffgeopy as ct
from IPython import embed


# enable latex printing in ipython qtconsole
qt = 1
if qt:
    from sympy import init_printing
    from IPython.display import display
    init_printing()

# specify how to optimize right pseudo inverse
#pinv_optimization = "count_ops" # options: "free_symbols", "count_ops"
pinv_optimization = "free_symbols" # options: "free_symbols", "count_ops"
#pinv_optimization = "none"


# read system from file: -----------------------------------------------

from system1 import *




string_line = "##################################################################\n\n"

def left_ortho_complement(matrix):
    assert has_left_ortho_complement(matrix), "There is no leftorthocomplement!"

    v = sp.simplify( st.nullspaceMatrix(matrix.T).T )
    assert is_zero_matrix(v*matrix), "Leftorthocomplement is not correct."
    return v

def has_left_ortho_complement(matrix):
    r = matrix.rank()
    m, n = matrix.shape
    if(n>=m):
        return True if (r<m) else False
    return True
    
def right_ortho_complement(matrix):
    assert has_right_ortho_complement(matrix), "There is no rightorthocomplement!"

    v = sp.simplify( st.nullspaceMatrix(matrix) )
    assert is_zero_matrix(matrix*v), "Leftorthocomplement is not correct."
    return v

def has_right_ortho_complement(matrix):
    r = matrix.rank()
    m, n = matrix.shape
    if(m>=n):
        return True if (r<n) else False
    return True

def left_pseudo_inverse(matrix):
    m, n = matrix.shape
    r = matrix.rank()
    assert r==n, "Matrix does not have full column rank!"

    lpinv = sp.simplify( (matrix.T * matrix).inv() * matrix.T )
    assert is_unit_matrix(lpinv*matrix), "Leftpseudoinverse is not correct."
    return lpinv

def get_column_index_of_matrix(matrix, column):
    m1, n1 = matrix.shape
    for index in xrange(n1):
        if matrix.col(index)==column:
            return index

def remove_zero_columns(matrix):
    m, n = matrix.shape
    M = sp.Matrix([])

    for i in xrange(n):
        if not is_zero_matrix(matrix.col(i)):
            M = st.concat_cols(M, matrix.col(i))
    return M

def is_linearly_independent(matrix, column_vector):
    # TODO: mit zahleneinsatz testen! st.subs_random_numbers()
    # dann mit rank und für beispiele ein vernünftiges epsilon finden
    m, n = matrix.shape
    for i in xrange(n):
        tmp = st.concat_cols(matrix.col(i),column_vector)
        if tmp.rank()!=2:
            return False
    return True

def matrix_to_vectorlist(matrix):
    m, n = matrix.shape
    vectors = []
    for i in xrange(n):
        vectors.append(matrix.col(i))
    return vectors

def count_zero_entries(vector):
    """ returns an integer with number of zeros
    """
    number = 0
    # check if zeros exist, if so: count them
    atoms = vector.atoms()
    if 0 in atoms:
        for i in xrange(len(vector)):
            if vector[i]==0:
                number+=1
    return number

def nr_of_ops(vector):
    m = vector.shape[0]
    ops = 0
    for i in xrange(m):
        ops += vector[i].count_ops()
    return ops

def right_pseudo_inverse(P):
    m0, n0 = P.shape

    # pick m0 of the simplest (=>(m0-1) or less entries are zero) lin.
    # independent columns of P:
    P_new = remove_zero_columns(P) # this step should not be necessary
    m1, n1 = P_new.shape

    list_of_cols = matrix_to_vectorlist(P_new)

    # sort by "complexity"
    if pinv_optimization=="free_symbols":
        cols_sorted_by_atoms = sorted( list_of_cols, key=lambda x: x.free_symbols, reverse=False)
    elif pinv_optimization=="count_ops":
        cols_sorted_by_atoms = sorted( list_of_cols, key=lambda x: nr_of_ops(x), reverse=False)
    elif pinv_optimization=="none":
        cols_sorted_by_atoms = list_of_cols

    # sort by number of zero entries
    colvecs_sorted = sorted( cols_sorted_by_atoms, key=lambda x: count_zero_entries(x), reverse=True)


    # pick m suitable column vectors and add to new matrix M
    M = colvecs_sorted[0]
    for j in xrange(len(colvecs_sorted)-1):
        column = colvecs_sorted[j+1]
        if is_linearly_independent(M,column):
            M = st.concat_cols(M,column)

        if M.rank()==m1:
            break

    # calculate R transformation matrix:
    R_tilde = sp.Matrix([])
    used_cols = []
    R = sp.Matrix([])
    for k in xrange(m1):
        new_column_index = k
        old_column_index = get_column_index_of_matrix(P,M.col(k))
        used_cols.append(old_column_index)
        tmp = sp.zeros(n0,1)
        tmp[old_column_index] = 1

        R_tilde = st.concat_cols(R_tilde,tmp)
    R=R_tilde

    # remainder columns of R matrix
    m2,n2 = R_tilde.shape
    for l in xrange(n0):
        if l not in used_cols:
            R_col = sp.zeros(n0,1)
            R_col[l] = 1
            R = st.concat_cols(R,R_col)

    m2, n2 = M.shape
    r2 = M.rank()
    assert m2==r2, "A problem occured in reshaping the Matrix"

    #M_inv = M.inv() #[m x m]
    M_det = M.berkowitz_det()
    M_inv = M.adjugate()/M_det
    #M.inv() #[m x m]
    p = n0-m1
    zero = sp.zeros(p,m1)
    Q = sp.simplify(M_inv)
    P_pinv = sp.simplify( R*st.concat_rows(Q,zero) )
    
    assert is_unit_matrix(sp.simplify(P*P_pinv)), "Rightpseudoinverse is not correct."
    return P_pinv


def diagonalize(matrix):
    """ TODO: ARE THERE ANY CASES WHERE THIS DOES NOT WORK???

        If r = matrix.rank(), then matrix can be diagonalized to
                        (I_r 0)
        matrix_diag =   (0   0)
    """
    m, n = matrix.shape
    r = matrix.rank()
    if (m>n) and (r!=n):
        raise NotImplementedError
    elif (n>m) and (r!=m):
        raise NotImplementedError
    
    
    if has_right_ortho_complement(matrix):
        roc = right_ortho_complement(matrix)
        rpinv = right_pseudo_inverse(matrix)
        transform_r = st.concat_cols(rpinv, roc)
        transform_l = sp.eye(matrix.shape[0])
        assert is_regular_matrix(transform_r), "Transformation matrix is not regular."
        diag = sp.simplify(matrix*transform_r)

    elif has_left_ortho_complement(matrix):
        loc = left_ortho_complement(matrix)
        lpinv = left_pseudo_inverse(matrix)
        transform_l = st.concat_rows(lpinv, loc)
        transform_r = sp.eye(matrix.shape[1])
        assert is_regular_matrix(transform_l), "Transformation matrix is not regular."
        diag = sp.simplify(transform_l*matrix)
    
    else:
        # quadratisch, voller rang
        transform_r = matrix.inv()
        transform_l = sp.eye(m)
        diag = sp.simplify(matrix*transform_r)
    
    return transform_l, diag, transform_r

def end_condition(Bi):
    """ The algorithm ends if Bi has full row rank.
    """
    n, p = Bi.shape
    return True if (n == Bi.rank()) else False

def is_zero_matrix(matrix):
    matrix = sp.simplify(matrix)
    m, n = matrix.shape
    zero_matrix = sp.zeros(m,n)
    return True if (matrix == zero_matrix) else False

def is_unit_matrix(matrix):
    matrix = sp.simplify(matrix)
    m, n = matrix.shape
    assert is_square_matrix(matrix), "Matrix is not a square matrix."
    unit_matrix = sp.eye(m)
    return True if (matrix == unit_matrix) else False

def is_square_matrix(matrix):
    m, n = matrix.shape
    return True if (m==n) else False

def is_regular_matrix(matrix):
    assert is_square_matrix(matrix), "Matrix is not a square matrix."
    m, n = matrix.shape
    r = matrix.rank()
    return True if (r == m) else False

def has_full_row_rank(matrix):
    m, n = matrix.shape
    r = matrix.rank()
    return True if (m==r) else False

def outlier(Bi):
    """ sonderfall 4.7
        tritt auf wenn Bi linear abhängige spalten aufweisst bzw.
        keinen vollen spaltenrang hat
    """
    n, p = Bi.shape
    return True if (Bi.rank() < p) else False

def assemble_dual_basis(Q):#, basis, basis_1forms):
    # einsform der dualen basis:
    p1, p2 = Q.shape
    l1 = len(basis)

    # w ist a vector, but basis_1form cannot currently be used in sympy
    # matrices, so a list will hold the entries of the column vector w:
    # create list with length p1 (p1 rows)
    w_T = [0] * p1
    for j in xrange(p1):
        # initialize with "empty" 1forms
        w_T[j] = 0 * basis_1form[0]
        for i in xrange(p2):
            w_T[j] = w_T[j] + Q.row(j)[i] * basis_1form[i]
    return w_T

def is_integrable_1(w):
    return all(x.d.is_zero() for x in w)

def is_integrable_2(w, i):
    for j in xrange(len(w)):
        if (w[i].d^w[j]).is_zero():
            print "w[" + str(i) + "].d^w[" + str(j) + "]=0"
            return True
        elif (w[i].d^w[j].d).is_zero() and j!=i:
            print "w[" + str(i) + "].d^w[" + str(j) + "].d=0"
            return True
    return False

def integrability_conditions(w):
    j = len(w)
    y = []
    for i in xrange(j):
        if (w[i].d).is_zero():
            print "w[" + str(i) + "].d = 0"
            print "<-> Integrabilitätsbedingung für w[" + str(i) + "] erfüllt."
            print "    Der Flache Ausgang wurde berechnet:"
            y=0
            for j in xrange(0,len(basis)):
                y += sp.integrate(w[i].coeff[j],basis[j])
            print "y" + str(i+1) + " ="
            print_nicely(y)
            print "\n\n"
        elif is_integrable_2(w,i):
            print "<-> Integrabilitätsbedingung für w[" + str(i) + "] erfüllt.\n\n"
        else:
            print "Integrabilitätsbedingung für w[" + str(i) + "] nicht erfüllt.\n\n"
    
    #if is_integrable_1(w):
    #    integrate_1form(w)

def integrate_1form(w):
    y = []
    for i in xrange(0,len(w)):
        y.append(0)
        for j in xrange(0,len(basis)):
            y[i] += sp.integrate(w[i].coeff[j],basis[j])

    print "Der flache Ausgang wurde berechnet:\n"
    print_equation_from_list(y, "y")

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

def print_matrix(name, i, supplement, matrix):
    print name + str(i) + supplement + " [" + str(matrix.shape[0]) + " x "  + str(matrix.shape[1]) + "] = "
    print_nicely(matrix)

def print_nicely(formula):
    if qt:
        display(formula)
    else:
        sp.pprint(formula)
    print("\n")

def print_line():
    print string_line

def print_equation_from_list(f, var, sympy=True):
    if sympy==False:
        for k in xrange(0, len(f)):
            print str(var) + str(k+1) + " = " + str(f[k]) + "\n"
    else:
        for k in xrange(0, len(f)):
            print str(var) + str(k+1) + " ="
            print_nicely(f[k])    

def check_row_compatibility(*args):
    """ checks if all matrices have same row dimension
    """
    return all(x.shape[0] == args[0].shape[0] for x in args)

def check_col_compatibilty(*args):
    """ checks if all matrices have same column dimension
    """
    return all(x.shape[1] == args[0].shape[1] for x in args)

def numerically_is_unit_matrix(matrix):
    m, n = matrix.shape
    r = matrix.rank()
    if (m!=n) or (m!=r):
        return False
    

def Zi_left_pinv_with_restrictions(Zi, P1i_tilde_right_ortho, P1i_right_pseudo_inv):
    """ Given a matrix Zi, this function calculates a matrix such that:
            Zi_left_pinv * Zi                    = I
            Zi_left_pinv * P1i_tilde_right_ortho = 0
            Zi_left_pinv * P1i_right_pseudo_inv  = 0
    """
    assert check_row_compatibility(Zi, P1i_tilde_right_ortho, P1i_right_pseudo_inv),\
        "Matrices do not match in row dimension."

    C = st.concat_cols(Zi, P1i_tilde_right_ortho, P1i_right_pseudo_inv)

    assert is_regular_matrix(C), "C is not a regular matrix"
    C_det = C.berkowitz_det()
    C_inv = C.adjugate()/C_det
    C_inv = sp.simplify(C.inv())

    m, n = Zi.shape
    Zi_left_pinv = sp.Matrix([])
    for i in xrange(0,n):
        Zi_left_pinv = st.concat_rows(Zi_left_pinv,C_inv.row(i))

    o, p = Zi_left_pinv.shape
    assert o==n and p==m, "There must have been a problem with the\
                            computation of Zi_left_pinv"


    # TODO: use st.subs_random_numbers() and np.allclose() instead!!!
    # muss alles erfüllt sein, trotzdem sind einige dieser ausdrücke
    # zu kompliziert als dass sp.simplify sie zur einheitsmatrix bzw.
    # nullmatrix verinfachen könnte. deshalb wird (vorerst) auf die
    # folgenden asserts verzichtet:

    #assert st.subs_random_numbers(Zi_left_pinv*Zi), "Zi_left_pinv is wrong"

    #assert is_unit_matrix(Zi_left_pinv*Zi), "Zi_left_pinv is wrong"
    #assert is_zero_matrix(Zi_left_pinv*P1i_tilde_right_ortho), "Zi_left_pinv is wrong"
    #assert is_zero_matrix(Zi_left_pinv*P1i_right_pseudo_inv), "Zi_left_pinv is wrong"
    return Zi_left_pinv


def first_substitution(P2i, P1i, P0i, i):

    P2i_roc = right_ortho_complement(P2i)
    P2i_rpinv = right_pseudo_inverse(P2i)

    P2i_dot = st.perform_time_derivative(P2i, diffvec_x)
    P2i_roc_dot = st.perform_time_derivative(P2i_roc, diffvec_x)
    P2i_roc_ddot = st.perform_time_derivative(P2i_roc, diffvec_x, order=2)

    Ai = sp.simplify( (P1i - 2*P2i_dot)*P2i_roc )
    Bi = sp.simplify( P2i*P2i_roc_ddot + P1i*P2i_roc_dot + P0i*P2i_roc )
    #Bi = sp.simplify( (P0i - P1i_dot)*P1i_roc )

    # print matrices
    print_matrix("P2",i,"_roc",P2i_roc)
    print_matrix("P2",i,"_rpinv",P2i_rpinv)
    print_matrix("P2",i,"_dot",P2i_dot)
    print_matrix("P2",i,"_roc_dot",P2i_roc_dot)
    print_matrix("P2",i,"_roc_ddot",P2i_roc_ddot)

    print_matrix("A",i,"",Ai)
    print_matrix("B",i,"",Bi)

    return Ai, Bi, P2i_rpinv,P2i_roc

def second_substitution(Ai, Bi, i):
    Ai_roc = right_ortho_complement(Ai)
    Ai_dot = st.perform_time_derivative(Ai, diffvec_x)
    #Ai_roc_dot = st.perform_time_derivative(Ai_roc, diffvec_x)
    
    Ci = sp.simplify( (Bi - Ai_dot)*Ai_roc )

    # print matrices
    print_matrix("A",i,"_roc",Ai_roc)
    print_matrix("A",i,"_dot",Ai_dot)
    print_matrix("C",i,"",Ci)
    return Ci
    


def fourseven(Ai, Bi, P1i_roc, P1i_rpinv, i):
    K2 = right_ortho_complement(Bi)
    print_nicely(K2)
    #K1 = right_ortho_complement(K2.T)
    K1 = right_pseudo_inverse(Bi)
    print_nicely(K1)

    K = st.concat_cols(K1, K2)

    assert is_regular_matrix(K), "K is not a regular matrix."
    
    #embed()

    Bi_tilde = sp.simplify(Bi*K1) # unit matrix
    P1i_tilde_roc = sp.simplify( P1i_roc*K1 )
    Zi = sp.simplify( P1i_roc*K2 )

    Zi_lpinv = Zi_left_pinv_with_restrictions(Zi, P1i_tilde_roc, P1i_rpinv)
    
    # print matrices
    print_matrix("K",i,"",K)
    print_matrix("B",i,"_tilde",Bi_tilde)
    print_matrix("Z",i,"",Zi)
    print_matrix("Z",i,"_lpinv",Zi_lpinv)

    return Ai, Bi_tilde, Zi_lpinv


def elimination(Ai, Bi, i):
    Bi_loc = left_ortho_complement(Bi)
    P1i_new = Bi_loc
    P0i_new = Bi_loc*Ai

    # print matrices
    print_matrix("B",i,"_loc",Bi_loc)

    return P1i_new, P0i_new


def calculate_Q_matrix():
    # TODO: unit test!!!
    
    Q = P1i_stack[0]
    # stack to hold P1i multiplications [P10,P11*P10,P12*P11P10,...]
    multiplication_stack = [Q]
    P1i_keys = P1i_stack.keys()
    for i in xrange(1, len(P1i_keys)):
        key = P1i_keys[i]
        P1i = P1i_stack[key]
        Q = P1i*Q
        multiplication_stack.append(Q)


    Zi_lpinv_stack_keys = Zi_lpinv_stack.keys()
    for j in xrange(len(Zi_lpinv_stack.keys())):
        key = Zi_lpinv_stack_keys[j]
        if key==0:
            Q_tilde = Zi_lpinv_stack[key]
        else:
            P = multiplication_stack[key-1] # P1i*P1(i-1)*...*P10 (wobei i=key)
            Q_tilde = Zi_lpinv_stack[key]*P

        Q = st.concat_rows(Q, Q_tilde)

    print "Q [" + str(Q.shape[0]) + " x " + str(Q.shape[1]) + "] ="
    print_nicely(Q)   
    return Q

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################


vec_xdot = st.perform_time_derivative(vec_x,vec_x)
vec_xddot = st.perform_time_derivative(vec_xdot,vec_xdot)
diffvec_x = st.concat_rows(vec_x, vec_xdot, vec_xddot)


print "x ="; print_nicely(vec_x)
print "\n\n xdot ="; print_nicely(vec_xdot)


# exterior derivative of F_eq:
# TODO: Not elegant, could be combined!
try:
    P2i # in case the system is given by the matrices P1i and P0i
except NameError:
    print "\n\n0 = F(x,xdot) ="; print_nicely(F_eq)
    P2i = sp.Matrix([])
    for i in xrange(len(vec_xdot)):
        vector = sp.simplify( F_eq.diff(vec_xdot[i]) )
        P2i = st.concat_cols(P2i, vector)

    P1i = sp.Matrix([])
    for i in xrange(len(vec_xdot)):
        vector = sp.simplify( F_eq.diff(vec_xdot[i]) )
        P1i = st.concat_cols(P1i, vector)

    P0i = sp.Matrix([])
    for i in xrange(len(vec_x)):
        vector = sp.simplify( F_eq.diff(vec_x[i]) )
        P0i = st.concat_cols(P0i, vector)
print "\n\n\n"

########################################################################
########################################################################
########################################################################
########################################################################

P2i_stack = {}
Zi_lpinv_stack = {}

i = 0
while 1:
    print "i = " + str(i) + " " + string_line[len(str(i))+5:]
    print_matrix("P2",i,"",P2i)
    print_matrix("P1",i,"",P1i)
    print_matrix("P0",i,"",P0i)
    # kann nur am anfang passieren:
    assert has_full_row_rank(P2i), "P20 does not have full row rank. There \
                                    must be algebraic equations."

    # add to stack
    P2i_stack[i] = P2i

    # 1. reduktionsschritt
    #Ai, Bi, P1i_roc, P1i_rpinv = first_substitution(P2i, P1i, P0i, i)
    Ai, Bi, P2i_rpinv, P2i_roc = first_substitution(P2i, P1i, P0i, i)

    K1 = P2i_rpinv
    K2 = P2i_roc
    K = st.concat_cols(K1, K2)
    K_inv = K.inv()

    embed()
    Ci = second_substitution(Ai, Bi, i)
    embed()

    assert not is_zero_matrix(Bi), "System ist not flat!"

    if outlier(Bi):
        print "Sonderfall"
        Ai, Bi, Zi_lpinv = fourseven(Ai, Bi, P1i_roc, P1i_rpinv, i)

        # add to stack
        Zi_lpinv_stack[i] = Zi_lpinv

    if end_condition(Bi):
        print_line()
        print("Algorithmus am Ende")
        break

    # 2. eliminationsschritt
    P1i, P0i = elimination(Ai, Bi, i)

    i += 1

# berechne neues Q
Q = calculate_Q_matrix()

# TODO: not the most elegant way?
# check highest order of derivatives
highest_order = 0
for n in Q.atoms():
    if hasattr(n, "difforder"):
        if n.difforder>highest_order:
            highest_order = n.difforder

# generate vector with vec_x and its derivatives up to highest_order
new_vector = vec_x
for index in xrange(1, highest_order+1):
    vec_x_ndot = st.perform_time_derivative(vec_x, vec_x, order=index)
    new_vector = st.concat_rows( new_vector, vec_x_ndot )

# generate basis_1form up to this order
basis, basis_1form = ct.diffgeo_setup(new_vector)

# calculate w
w = assemble_dual_basis(Q)#, basis, basis_1form)

print_equation_from_list(w, "w", sympy=False)
print_line()

#try_integrability_conditions(w)
integrability_conditions(w)

print_line()



