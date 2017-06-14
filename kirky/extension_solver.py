from cvxopt import matrix, solvers
from fractions import Fraction
from numpy.linalg import matrix_rank, qr
from numpy import array
from copy import deepcopy 
"""
For information on the extension being used see: http://cvxopt.org/userguide/coneprog.html#linear-programming
"""
def resolve_repeating_fraction(rf):
    ACCUR = 13
    for q in range(1, 14):
        for k in range(1, q + 1):
            attempt = Fraction(rf).limit_denominator(10**q*10**(k-1))
            print (float(attempt) -rf)
            if abs(float(attempt) - rf) < 10**-12:
                return attempt


def float_matrix(A):
    F = []
    for row in A:
        F.append([float(e) for e in row])
    return F


def reduce_matrix(A):
    new_A = [A[0]]
    for i in range(1, len(A)):
        if matrix_rank(array(float_matrix(new_A) + float_matrix([A[i]]))) > matrix_rank(array(float_matrix(new_A))):
            new_A.append(A[i])
    return new_A        


def create_augmented_matrix(E):
    """
    This creates from matrix E a matrix with the sum condition 
    row.
    """
    A = deepcopy(E)
    m = len(A[0])
    n = len(A)
    for i in range(n):
        A[i].append(Fraction(0))
    sum_condition_row = [Fraction(1)] * m + [Fraction(-1)]
    A.append(sum_condition_row)
    return A


def get_b(Au):
    """
    This expects the augmented matrix not just plain old E
    
    returns a column vector of fractions
    """
    n = len(Au)
    b = [Fraction(0)] * (n-1) + [Fraction(1)]
    return b


def retrieve_components(Au):
    """
    This expects the augmented matrix, not just plain old E
    
    So cvxopt matrix expects columns. Rather every deepest list within a matrix constructor is
    a column. So we're going to need to extract the columns of E.
    """
    f_Au = float_matrix(Au)
    n = len(f_Au)
    m = len(f_Au[0])
    columns = [[f_Au[i][j] for i in range(n)] for j in range(m)]
    # next we add the identity bit at the end here to make sure we have a basis
    auxiliary_columns = []
    for i in range(n):
        new_column = [0.0] * n
        new_column[i] = 1.0
        auxiliary_columns.append(new_column)
    columns.extend(auxiliary_columns)
    # now we have our matrix A
    A = matrix(columns)
    # next we create our objective vector    
    c = matrix([0.0]*m + [1.0]*n)
    # next we create b
    b = [float(el) for el in get_b(Au)]
    b = matrix(b)
    # now we create G which is just going to be a negative identity
    columns_of_G = []
    for i in range(m+n):
        new_column = [0.0] * (n+m)
        new_column[i] = -1.0
        columns_of_G.append(new_column)
    G = matrix(columns_of_G)
    # and finally we have h
    h = matrix([0.0] * (n+m))
    # and we return them all
    return c, G, h, A, b


def get_solution(Au):
    c, G, h, A, b = retrieve_components(Au)
    dims = {
        'l': G.size[0],
        'q': [],
        's': []
    }
    solution = solvers.conelp(c, G, h, dims, A, b)
    return solution


def grab_exact_solution(A, basis):
    # this just pivots on the basis to find an exact solution
    rows_without_pivot = list(range(len(A)))
    c_to_r = {}
    for c in basis:
        found_a_pivot_row = False
        for x in range(len(rows_without_pivot)):
            if A[rows_without_pivot[x]][c] != 0:
                found_a_pivot_row = True
		break
	if found_a_pivot_row:
            r = rows_without_pivot.pop(x)
            c_to_r[c] = r
        else:
            print 'OMG!'
            #rint [A[i][c] for i in rows_without_pivot]
            continue
        pivot_el = A[r][c]
        A[r] = [el / pivot_el for el in A[r]]
        for i in range(len(A)):
            if i == r:
		continue
            coefficient = A[i][c]
            A[i] = [A[i][j] - coefficient * A[r][j] for j in range(len(A[i]))]
    # so at this point we've row reduced on the pivots so now we just 
    # need to get our solution
    solution = [Fraction(0.0)] * (len(A[0]) - 1) # the minus one is because we expect this to be augmented
    for c in c_to_r:
        r = c_to_r[c]
        solution[c] = A[r][-1]
    return solution

def rank_of_basis(A, basis):
    B = []
    for r in range(len(A)):
        B.append([float(A[r][c]) for c in basis])
    return matrix_rank(array(B))
	    

def check_solution(Au, ibos):
    # ibos - important_bits_of_solution
    print "INEXACT RESULTS OF SOLUTION"
    for i in range(len(Au)):
        print sum([float(Au[i][j]) * ibos[j] for j in range(len(ibos))])


def solve_kirky_with_extension(E):
    solution = grab_solution(E)
    if not solution:
        return
    else:
        return solution

def clean_solution(solution):
    def drop(e):
        if abs(e) < 10 ** -8:
            return 0.0
        else:
            return e
    return [drop(e) for e in solution]

def grab_solution(E):
    E = reduce_matrix(E)
    Au = create_augmented_matrix(E)
    solution = get_solution(Au)
    print solution['primal objective']
    if solution['primal objective'] > 10 ** -5: #TODO better decision making here
        return
    else:
        m = len(Au[0])
        important_bits_of_solution = [x for x in solution['x']][:m-1]
        c_solution = clean_solution(important_bits_of_solution)
        print c_solution
        return c_solution

"""
So it's definitely finding a good solution, and it's finding the right basis...
so either my exact solver is wrong, or there is something wrong in my logic. Let's
make a really damn simple exact solver then. 
"""
