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
import numpy as np
import symbtools as st

from IPython import embed as IPS

# specify how to optimize calculation of pseudo inverse / ortho complement
pinv_optimization = "free_symbols" # options: "free_symbols", "count_ops", "none"

roc = "alt" # options: "alt", "conv"
debug = False


def not_simplify(expr, **kwargs):
    return expr

custom_simplify = sp.simplify
custom_simplify = not_simplify


from util.timeout import get_timed_simplify

custom_simplify = get_timed_simplify(20)

st.nullspace_simplify_func = custom_simplify

# use prime number in subs_random_numbers()?
srn_prime = True

def left_ortho_complement(matrix):
    
    # commented out for performance
    #assert has_left_ortho_complement(matrix), "There is no leftorthocomplement!"

    v = custom_simplify( st.nullspaceMatrix(matrix.T).T )
    assert is_zero_matrix(v*matrix), "Leftorthocomplement is not correct."
    return v

def has_left_ortho_complement(matrix):
    r = st.rnd_number_rank(matrix)
    m, n = matrix.shape
    if(n>=m):
        return True if (r<m) else False
    return True
    
def right_ortho_complement(matrix):
    assert has_right_ortho_complement(matrix), "There is no rightorthocomplement!"

    if roc=="conv":
        v = custom_simplify( st.nullspaceMatrix(matrix) )
        assert is_zero_matrix(matrix*v), "Leftorthocomplement is not correct."
        return v
    else:
        return alternative_right_ortho_complement(matrix)

def has_right_ortho_complement(matrix):
    r = st.rnd_number_rank(matrix)
    m, n = matrix.shape
    if(m>=n):
        return True if (r<n) else False
    return True

def left_pseudo_inverse(matrix):
    m, n = matrix.shape
    r = st.rnd_number_rank(matrix)
    #assert r==n, "Matrix does not have full column rank!"
    assert not is_zero_matrix(matrix), "Matrix is zero matrix!"

    transposed_rpinv = right_pseudo_inverse(matrix.T)
    matrix_lpinv = transposed_rpinv.T

    assert is_unit_matrix(matrix_lpinv*matrix), "Leftpseudoinverse is not correct."

    return matrix_lpinv

def get_column_index_of_matrix(matrix, column):
    """ if V is a column vector of a matrix A, than this function returns
        the column index.
    """
    m1, n1 = matrix.shape
    for index in range(n1):
        if matrix.col(index)==column:
            return index
    return None

def remove_zero_columns(matrix):
    """ this function removes zero columns of a matrix
    """
    m, n = matrix.shape
    M = sp.Matrix([])

    for i in range(n):
        if not is_zero_matrix(matrix.col(i)):
            M = st.concat_cols(M, matrix.col(i))
    return M

def is_linearly_independent(matrix, column_vector):
    m, n = matrix.shape

    rank1 = st.rnd_number_rank(matrix)
    tmp = st.concat_cols(matrix, column_vector)
    rank2 = st.rnd_number_rank(tmp)
    
    assert rank2 >= rank1
        
    return rank2 > rank1

def matrix_to_vectorlist(matrix):
    m, n = matrix.shape
    vectors = []
    for i in range(n):
        vectors.append(matrix.col(i))
    return vectors

def count_zero_entries(vector):
    """ returns an integer with number of zeros
    """
    number = 0
    # check if zeros exist, if so: count them
    atoms = vector.atoms()
    if 0 in atoms:
        for i in range(len(vector)):
            if vector[i]==0:
                number+=1
    return number

def nr_of_ops(vector):
    m = vector.shape[0]
    ops = 0
    for i in range(m):
        ops += vector[i].count_ops()
    return ops

def regular_completion(matrix):
    m,n = matrix.shape
    r = matrix.rank()
    
    #~ IPS()
    assert m!=n, "There is no regular completion of a square matrix."

    if m<n:
        assert r==m, "Matrix does not have full row rank."
        A, B, V_pi = reshape_matrix_columns(matrix)
        zeros = sp.zeros(n-m,m)
        ones = sp.eye(n-m)
        S = st.col_stack(zeros,ones)
        completion = S*V_pi.T
        
        regular_matrix = st.row_stack(matrix,completion)
        assert st.rnd_number_rank(regular_matrix)==n, "Regular completion seems to be wrong."

    elif n<m:
        assert r==n, "Matrix does not have full column rank."
        A, B, V_pi = reshape_matrix_columns(matrix.T)
        zeros = sp.zeros(m-n,n)
        ones = sp.eye(m-n)
        S = st.col_stack(zeros,ones)
        completion = V_pi*S.T
        
        regular_matrix = st.col_stack(completion,matrix)
        assert st.rnd_number_rank(regular_matrix)==m, "Regular completion seems to be wrong."

    return completion

def reshape_matrix_columns(P):
    m0, n0 = P.shape

    # pick m0 of the simplest lin. independent columns of P:
    P_new = remove_zero_columns(P)
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
    for j in range(len(colvecs_sorted)-1):
        column = colvecs_sorted[j+1]
        if is_linearly_independent(A,column):
            A = st.concat_cols(A,column)

        if st.rnd_number_rank(A)==m1:
            break

    assert A.is_square
            
    # calculate transformation matrix R: -------------------------------
    #R_tilde = sp.Matrix([])
    used_cols = []
    R = sp.Matrix([])
    for k in range(m1):
        new_column_index = k
        old_column_index = get_column_index_of_matrix(P,A.col(k))
        used_cols.append(old_column_index)
        tmp = sp.zeros(n0,1)
        tmp[old_column_index] = 1

        #R_tilde = st.concat_cols(R_tilde,tmp)
        R = st.concat_cols(R,tmp)
    #R=R_tilde

    # remainder columns of R matrix
    #m2,n2 = R_tilde.shape
    m2,n2 = R.shape
    for l in range(n0):
        if l not in used_cols:
            R_col = sp.zeros(n0,1)
            R_col[l] = 1
            R = st.concat_cols(R,R_col)

    m3, n3 = A.shape
    r2 = st.rnd_number_rank(A)
    assert m3==r2, "A problem occured in reshaping the matrix."

    # calculate B matrix: ----------------------------------------------
    B = sp.Matrix([])
    tmp = P*R
    for i in range(n0-m0):
        B = st.concat_cols(B,tmp.col(m0+i))

    
    assert is_zero_matrix( (P*R) - st.concat_cols(A,B)), "A problem occured in reshaping the matrix."

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
    Q = custom_simplify(A_inv)
    P_pinv = custom_simplify( R*st.concat_rows(Q,zero) )

    assert is_unit_matrix(P*P_pinv), "Rightpseudoinverse is not correct."
    return P_pinv

def alternative_right_ortho_complement(P):
    m0, n0 = P.shape
    A, B, R = reshape_matrix_columns(P)

    #apparently this is quicker than A_inv = A.inv() #[m x m]:
    A_det = A.berkowitz_det()
    A_inv = A.adjugate()/A_det

    minusAinvB = custom_simplify(-A_inv*B)

    p = n0-m0
    unit_matrix = sp.eye(p)
    Q = custom_simplify(A_inv)
    P_roc = custom_simplify( R*st.concat_rows(minusAinvB, unit_matrix) )

    assert is_zero_matrix(P*P_roc), "Right orthocomplement is not correct."
    return P_roc


def Zi_left_pinv_with_restrictions(P1i_rpinv, P1i_tilde_roc, Zi):
    """ Given a matrix Zi, this function calculates a matrix such that:
            Zi_left_pinv * Zi            = I
            Zi_left_pinv * P1i_tilde_roc = 0
            Zi_left_pinv * P1i_rpinv     = 0
    """
    assert check_row_compatibility(P1i_rpinv, P1i_tilde_roc, Zi),\
        "Matrices do not match in row dimension."

    C = st.concat_cols(P1i_rpinv, P1i_tilde_roc, Zi)

    assert is_regular_matrix(C), "C is not a regular matrix"
    C_det = C.berkowitz_det()
    C_inv = C.adjugate()/C_det
    C_inv = custom_simplify(C_inv)

    m, n = Zi.shape
    Zi_left_pinv = sp.Matrix([])
    for i in range(m-n,m):
        Zi_left_pinv = st.concat_rows(Zi_left_pinv,C_inv.row(i))

    o, p = Zi_left_pinv.shape
    assert o==n and p==m, "There must have been a problem with the\
                            computation of Zi_left_pinv"

    assert is_unit_matrix(Zi_left_pinv*Zi), "Zi_left_pinv is wrong"
    assert is_zero_matrix(Zi_left_pinv*P1i_tilde_roc), "Zi_left_pinv is wrong"
    assert is_zero_matrix(Zi_left_pinv*P1i_rpinv), "Zi_left_pinv is wrong"

    return Zi_left_pinv


def is_symbolically_zero_matrix(matrix):
    matrix = custom_simplify(matrix)
    m, n = matrix.shape
    zero_matrix = sp.zeros(m,n)
    return True if (matrix == zero_matrix) else False

def is_zero_matrix(matrix):
    m_rand = st.subs_random_numbers(matrix, prime=srn_prime)

    for i in range(len(m_rand)):
        if not np.allclose(float(m_rand[i]), 0):
            return False
    return True

def is_symbolically_unit_matrix(matrix):
    matrix = custom_simplify(matrix)
    m, n = matrix.shape
    assert is_square_matrix(matrix), "Matrix is not a square matrix."
    unit_matrix = sp.eye(m)
    return True if (matrix == unit_matrix) else False

def is_unit_matrix(matrix):
    assert is_square_matrix(matrix), "Matrix is not a square matrix."
    m, n = matrix.shape
    m_rand = st.subs_random_numbers(matrix, prime=srn_prime)

    for i in range(len(m_rand)):
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
    r = st.rnd_number_rank(matrix)
    return True if (r == m) else False

def has_full_row_rank(matrix):
    m, n = matrix.shape
    r = st.rnd_number_rank(matrix)
    return True if (m==r) else False

def check_row_compatibility(*args):
    """ checks if all matrices have same row dimension
    """
    return all(x.shape[0] == args[0].shape[0] for x in args)

def check_col_compatibilty(*args):
    """ checks if all matrices have same column dimension
    """
    return all(x.shape[1] == args[0].shape[1] for x in args)
