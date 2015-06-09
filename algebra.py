# -*- coding: utf-8 -*-

import sympy as sp
import numpy as np
import symb_tools as st

# specify how to optimize calculation of pseudo inverse / ortho complement
pinv_optimization = "free_symbols" # options: "free_symbols", "count_ops", "none"

roc = "alt" # options: "alt", "conv"

def srank(matrix):
    """ Computes the rank for symbolic matrices.
    """
    m, n = matrix.shape
    matrix_rand = st.subs_random_numbers(matrix)
    r = matrix_rand.rank()

    if r == m:
        return r
    else:
        r_symb = matrix.rank()
        if r == r_symb:
            return r
        else:
            # TODO: was hier? bisher noch nicht aufgetreten. tests?
            pass

def left_ortho_complement(matrix):
    assert has_left_ortho_complement(matrix), "There is no leftorthocomplement!"

    v = sp.simplify( st.nullspaceMatrix(matrix.T).T )
    assert is_zero_matrix(v*matrix), "Leftorthocomplement is not correct."
    return v

def has_left_ortho_complement(matrix):
    r = srank(matrix)
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
    r = srank(matrix)
    m, n = matrix.shape
    if(m>=n):
        return True if (r<n) else False
    return True

def left_pseudo_inverse(matrix):
    m, n = matrix.shape
    r = srank(matrix)
    assert r==n, "Matrix does not have full column rank!"

    transposed_rpinv = right_pseudo_inverse(matrix.T)
    matrix_lpinv = transposed_rpinv.T

    assert is_unit_matrix(matrix_lpinv*matrix), "Leftpseudoinverse is not correct."

    return matrix_lpinv

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
        if srank(tmp)!=2:
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

def reshape_matrix_columns(P):
    m0, n0 = P.shape

    # pick m0 of the simplest (=>(m0-1) or less entries are zero) lin.
    # independent columns of P:
    P_new = remove_zero_columns(P) # normally, this step should not be necessary
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

    # pick m suitable column vectors and add to new matrix A: ----------
    A = colvecs_sorted[0]
    for j in xrange(len(colvecs_sorted)-1):
        column = colvecs_sorted[j+1]
        if is_linearly_independent(A,column):
            A = st.concat_cols(A,column)

        if srank(A)==m1:
            break

    # calculate transformation matrix R: -------------------------------
    R_tilde = sp.Matrix([])
    used_cols = []
    R = sp.Matrix([])
    for k in xrange(m1):
        new_column_index = k
        old_column_index = get_column_index_of_matrix(P,A.col(k))
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

    m2, n2 = A.shape
    #r2 = A.rank()
    r2 = srank(A)
    assert m2==r2, "A problem occured in reshaping the matrix."

    # calculate B matrix: ----------------------------------------------
    B = sp.Matrix([])
    tmp = P*R
    for i in xrange(n0-m0):
        B = st.concat_cols(B,tmp.col(m0+i))
    return A, B, R

def right_pseudo_inverse(P):
    """ Calculates a right pseudo inverse with as many zero entries as
        possible. Given a [m x n] matrix P, the algorithm picks m
        linearly independent column vectors of P to form a regular
        block matrix A such that Q = P*R = (A B) with the permuation matrix
        P. The right pseudo inverse of Q is
            Q_rpinv = ( A^(-1) )
                      (   0    )
        and the right pseudo inverse of P
            P_rpinv = R*( A^(-1) )
                        (   0    ).
        There is a degree of freedom in choosing the column vectors for
        A. At the moment this can be done with respect to "count_ops"
        (minimizes number of operations) or "free_symbols" (minimizes the
        number of symbols).
    """
    m0, n0 = P.shape
    A, B, R = reshape_matrix_columns(P)

    #apparently this is quicker than A_inv = A.inv() #[m x m]:
    A_det = A.berkowitz_det()
    A_inv = A.adjugate()/A_det

    p = n0-m0
    zero = sp.zeros(p,m0)
    Q = sp.simplify(A_inv)
    P_pinv = sp.simplify( R*st.concat_rows(Q,zero) )

    assert is_unit_matrix(P*P_pinv), "Rightpseudoinverse is not correct."
    return P_pinv

def alternative_right_ortho_complement(P):
    m0, n0 = P.shape
    A, B, R = reshape_matrix_columns(P)

    #apparently this is quicker than A_inv = A.inv() #[m x m]:
    A_det = A.berkowitz_det()
    A_inv = A.adjugate()/A_det

    minusAinvB = sp.simplify(-A_inv*B)

    p = n0-m0
    unit_matrix = sp.eye(p)
    Q = sp.simplify(A_inv)
    P_roc = sp.simplify( R*st.concat_rows(minusAinvB, unit_matrix) )

    assert is_zero_matrix(P*P_roc), "Right orthocomplement is not correct."
    return P_roc


#def diagonalize(matrix):
    #""" TODO: ARE THERE ANY CASES WHERE THIS DOES NOT WORK???

        #If r = matrix.rank(), then matrix can be diagonalized to
                        #(I_r 0)
        #matrix_diag =   (0   0)
    #"""
    #m, n = matrix.shape
    #r = matrix.rank()
    #if (m>n) and (r!=n):
        #raise NotImplementedError
    #elif (n>m) and (r!=m):
        #raise NotImplementedError
    
    
    #if has_right_ortho_complement(matrix):
        #roc = right_ortho_complement(matrix)
        #rpinv = right_pseudo_inverse(matrix)
        #transform_r = st.concat_cols(rpinv, roc)
        #transform_l = sp.eye(matrix.shape[0])
        #assert is_regular_matrix(transform_r), "Transformation matrix is not regular."
        #diag = sp.simplify(matrix*transform_r)

    #elif has_left_ortho_complement(matrix):
        #loc = left_ortho_complement(matrix)
        #lpinv = left_pseudo_inverse(matrix)
        #transform_l = st.concat_rows(lpinv, loc)
        #transform_r = sp.eye(matrix.shape[1])
        #assert is_regular_matrix(transform_l), "Transformation matrix is not regular."
        #diag = sp.simplify(transform_l*matrix)
    
    #else:
        ## quadratisch, voller rang
        #transform_r = matrix.inv()
        #transform_l = sp.eye(m)
        #diag = sp.simplify(matrix*transform_r)
    
    #return transform_l, diag, transform_r


def is_symbolically_zero_matrix(matrix):
    matrix = sp.simplify(matrix)
    m, n = matrix.shape
    zero_matrix = sp.zeros(m,n)
    return True if (matrix == zero_matrix) else False

def is_zero_matrix(matrix):
    m_rand = st.subs_random_numbers(matrix)

    for i in xrange(len(m_rand)):
        if not np.allclose(float(m_rand[i]), 0):
            return False
    return True

def is_symbolically_unit_matrix(matrix):
    matrix = sp.simplify(matrix)
    m, n = matrix.shape
    assert is_square_matrix(matrix), "Matrix is not a square matrix."
    unit_matrix = sp.eye(m)
    return True if (matrix == unit_matrix) else False

def is_unit_matrix(matrix):
    assert is_square_matrix(matrix), "Matrix is not a square matrix."
    m, n = matrix.shape
    m_rand = st.subs_random_numbers(matrix)

    for i in xrange(len(m_rand)):
        if i%(m+1)==0:
            if not np.allclose(float(m_rand[i]), 1):
                return False
        else:
            if not np.allclose(float(m_rand[i]), 0):
                return False
    return True

def is_square_matrix(matrix):
    m, n = matrix.shape
    return True if (m==n) else False

def is_regular_matrix(matrix):
    assert is_square_matrix(matrix), "Matrix is not a square matrix."
    m, n = matrix.shape
    #r = matrix.rank()
    r = srank(matrix)
    return True if (r == m) else False

def has_full_row_rank(matrix):
    m, n = matrix.shape
    #r = matrix.rank()
    r = srank(matrix)
    return True if (m==r) else False

def check_row_compatibility(*args):
    """ checks if all matrices have same row dimension
    """
    return all(x.shape[0] == args[0].shape[0] for x in args)

def check_col_compatibilty(*args):
    """ checks if all matrices have same column dimension
    """
    return all(x.shape[1] == args[0].shape[1] for x in args)

