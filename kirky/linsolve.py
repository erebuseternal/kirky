from sklearn.svm import SVC
from fractions import Fraction
import numpy as np
from numpy.linalg import svd
from math import sqrt
from .splitter import Splitter
from sympy import Matrix, Rational

def reflect(_list):
    return [-a for a in _list]


def dot(vector1, vector2):
    return sum([vector1[i] * vector2[i] for i in range(len(vector1))])


def normalize(vector):
    length = sqrt(sum([element ** 2 for element in vector]))
    if length == 0.0:
        return [element for element in vector]
    else:
        return [element / length for element in vector]


def nullspace(matrix, atol=10**-13, rtol=0):
    """
    nullspace(matrix, atol=10**-14, rtol=0) - This function takes a numpy matrix and returns the nullspace
    matrix for the input. The code was copied from http://scipy-cookbook.readthedocs.io/items/RankNullspace.html
    Check that site out if you want to know more.

    But just in case the site disappears, this is what the site says:
        matrix should be at most 2-D.  A 1-D array with length k will be treated
            as a 2-D with shape (1, k)
        atol : float
            The absolute tolerance for a zero singular value.  Singular values
            smaller than `atol` are considered to be zero.
        rtol : float
            The relative tolerance.  Singular values less than rtol*smax are
            considered to be zero, where smax is the largest singular value.
    """
    matrix = np.atleast_2d(matrix)
    u, s, vh = svd(matrix)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns


def shave(matrix, tol=10**-12):
    """
    shave(matrix, tol=10**-13) - This function takes a matrix and goes through each of its elements
    using python's Fraction to find the closes number that does not have any part smaller than the
    input tolerance (tol). It does this by converting the element to a Fraction where the denominator
    is limited to being no greater than 1 / tol. It then spits out the resulting matrix.
    """
    rows = []
    max_denominator = int(1.0 / tol)
    for row in matrix:
        new_row = [float(Fraction(element).limit_denominator(max_denominator)) for element in row]
        rows.append(new_row)
    return np.array(rows)

def zero_shave(matrix):
    rows = []
    for row in matrix:
        new_row = []
        for e in row:
            if abs(e) <= 10 ** -12:
                new_row.append(0.0)
            else:
                new_row.append(e)
        rows.append(new_row)
    return np.array(rows)


def positive(matrix):
    """
    positive(matrix) - Given a numpy matrix, this function returns its positive if it exists and None otherwise
    """
    # we transpose the matrix so we can easily grab its columns
    trans = np.transpose(matrix)
    # then we grab the columns
    columns = [row for row in trans]
    # next we use the divider to try to divide the columns
    splitter = Splitter(columns)
    norm = splitter.split()
    print norm
    # we check to make sure the norm exists and thereby division was successful
    if norm is None:
        return None
    # if it was we then create a numpy row vector from the norm
    norm = np.array([[element for element in norm]])
    # and we dot it with the original matrix
    positive = np.dot(norm, matrix)
    # finally we shave the error from the positive and return it
    return shave(positive)


def positive_null_space_vector(matrix, verbose=True):
    """
    positive_nullspace_vector(matrix) - This function takes a matrix and tries to find
    a nullspace vector that has only positive entries. If it finds one it returns that
    vector, if it doesn't it returns None. (The vector is returned as a column vector)
    """
    # first we get the nullspace of the matrix
    N = nullspace(matrix)
    # we check to see if the nullspace is empty
    if N.shape[0] == 0:
        # no positive to be found here
        if verbose: print 'nullspace is empty'
        return None
    # we shave off any error
    N = zero_shave(N)
    # now to generate the positive we need to transpose N, because right not the nullspace vectors
    # are its columns and we need them to be rows for this to work properly
    N = np.transpose(N)
    # now we can try to get the positive for this matrix
    pos = positive(N)
    # we check to see if a positive was found
    if pos is not None:
        # we found it! but it's gonna be a row vectors and nullspace vectors should be columns
        # so we'll transpose it before we send it on its merry way
        pos = np.transpose(pos)
        if verbose: print 'positive found'
        return pos
    else:
        # in this case no positive was found
        if verbose: print 'no positive found'
        return None

def exact_positive_null_space_vector(matrix, verbose=True):
    N = matrix.nullspace()
    if len(N) == 0:
        print 'no nullspace found'
        return None
    rows = []
    # we grab the rows of this matrix as our vectors
    for i in range(N[0].shape[0]):
        row = []
        for column in N:
            row.append(Fraction(column[i,0].p, column[i,0].q))
        rows.append(row)
    splitter = Splitter(rows, exact=True)
    norm = splitter.split()
    if norm is None:
        print 'no norm found'
        return None
    # now we make the transpose of N
    rows = []
    for column in N:
        row = [e for e in column]
        rows.append(row)
    trans = Matrix(rows)
    symp_norm = Matrix([[Rational(e.numerator, e.denominator) for e in norm]])
    result = symp_norm * trans
    return result.T
