from cvxopt import matrix, solvers
from fractions import Fraction
from numpy.linalg import matrix_rank, qr
from numpy import array
from copy import deepcopy 
"""
For information on the extension being used see: http://cvxopt.org/userguide/coneprog.html#linear-programming
"""
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
    E = reduce_matrix(E)
    Au = create_augmented_matrix(E)
    solution = grab_solution(Au)
    print solution
    if not solution:
        return
    else:
        return solution[:-1]

def grab_solution(E):
    Au = create_augmented_matrix(E)
    print 'NEW NUMBER IF COLS: %s' % len(Au[0])
    solution = get_solution(Au)
    print solution['primal objective']
    if solution['primal objective'] > 10 ** -5: #TODO better decision making here
        return
    else:
        # now I want to grab the biggest n where n is the number of rows in Au...
        # so we're getting a basis w.r.t to our huge system... So that's what
        # we really need to get at.
        m = len(Au[0])
        print solution['x']
        important_bits_of_solution = [x for x in solution['x']][:m]
        check_solution(Au, important_bits_of_solution)
        n = len(Au)
        maxes = set()
        for i in range(n):
            _max = max(important_bits_of_solution)
            for j in range(len(important_bits_of_solution)):
                if important_bits_of_solution[j] == _max:
                    important_bits_of_solution.pop(j)
                    break
            maxes.add(_max)
        print 'MAXES: %s' % maxes
        basis = []
        for i in range(m):
            if solution['x'][i] in maxes:
                basis.append(i)
        # now we create an augmented, augmented matrix (we're adding b to the end of Au
	print 'RANK: %s' % matrix_rank(array(float_matrix(Au)))
        print 'N: %s' % len(Au)
        print 'M: %s' % m
        print 'IN BASIS: %s' % len(basis)
        if rank_of_basis(Au, basis) < matrix_rank(array(float_matrix(Au))):
            basis = [c for c in basis if c != (m - 1)]
            B = []
    	    for r in range(len(Au) - 1):
                B.append([Au[r][c] for c in basis])
            print '%s NUMBER OF COLS' % len(B[0])
            B = reduce_matrix(B)
            print 'recursing now'
            small_solution = grab_solution(B)
            print len(basis)
            print len(small_solution)
            solution = [Fraction(0)] * (m - 1)
            print len(solution)
            for x in range(len(basis)):
                print basis[x]
                solution[basis[x]] = small_solution[x]
            return solution
        b = get_b(Au)
        for i in range(n):
            Au[i].append(b[i])    # we augement each row with the corresponding element from b
        cleaned_solution = grab_exact_solution(Au, basis)[:-1]
        print cleaned_solution
        return cleaned_solution

"""
So it's definitely finding a good solution, and it's finding the right basis...
so either my exact solver is wrong, or there is something wrong in my logic. Let's
make a really damn simple exact solver then. 
"""
