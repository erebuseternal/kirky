from cvxopt import matrix, solvers
from fractions import Fraction
from numpy.linalg import matrix_rank, qr
from numpy import array
from copy import deepcopy 
"""
For information on the extension being used see: http://cvxopt.org/userguide/coneprog.html#linear-programming
"""
########################## test_reduce_matrix #########################
def float_matrix(A):
    F = []
    for row in A:
        F.append([float(e) for e in row])
    return F


def conseq_ending_zeros(row):
    conseq_zeros = 0
    for el in list(row)[::-1]:
        if abs(el) < 10**-12:
            conseq_zeros += 1
        else:
            break
    return conseq_zeros

"""
def reduce_matrix(A):
    A = float_matrix(A)
    # we want to get a set of dependent columns 
    # out of there. QR decomp lets us get rid of 
    # dependent columns so transpose will work
    A = array(A)
    T = A.T
    Q, R = qr(T)
    # now we just need to find the columns of R that 
    # give us an triangular matrix these indices will
    # be the indices of our linearly dependent rows
    RT = R.T
    taken = set(range(RT.shape[1]))
    good_indices = list()
    new_A = []
    for i in range(RT.shape[0]):
        row = RT[i]
        conseq_zeros = conseq_ending_zeros(row)
        if conseq_zeros == 3:
            print row
        if conseq_zeros in taken:
            taken.remove(conseq_zeros) 
            new_A.append([el for el in A[i]])
    ############ TODO ############
    bad_rows = []
    print R
    for i in range(RT.shape[0]):
        if abs(RT[i][i]) < 10 ** -12:
            bad_rows.append(i)  
    new_A = []
    for i in range(len(A)):
        if i not in bad_rows:
            new_A.append([el for el in A[i]])
    ##############################
    return new_A
"""
def reduce_matrix(A):
    A = float_matrix(A)
    new_A = [A[0]]
    for i in range(1, len(A)):
        if matrix_rank(array(new_A + [A[i]])) > matrix_rank(array(new_A)):
            new_A.append(A[i])
    return new_A        

def reduce_matrix_exact(A):
    new_A = [A[0]]
    for i in range(1, len(A)):
        if matrix_rank(array(float_matrix(new_A) + float_matrix([A[i]]))) > matrix_rank(array(float_matrix(new_A))):
            new_A.append(A[i])
    return new_A        

########################## test_create_augmented_matrix #########################
def create_augmented_matrix(E):
    A = deepcopy(E)
    m = len(A[1])
    n = len(A)
    for i in range(n):
        A[i].append(0.0)
    sum_condition_row = [1.0] * m + [-1.0]
    A.append(sum_condition_row)
    return A
    
def retrieve_components(E):
    """
    So cvxopt matrix expects columns. Rather every deepest list within a matrix constructor is
    a column. So we're going to need to extract the columns of E.

    The extra 1 at the end is for the sum condition row
    """
    # TODO handle matrices as columns in order to be more efficient
    num_rows = len(E)
    num_columns = len(E[0])
    columns = [[float(E[i][j]) for i in range(num_rows)] + [1.0] for j in range(num_columns)]
    """
    We add the column for the sum condition column
    """
    sum_condition_column = [0.0] * num_rows + [-1.0]
    columns.append(sum_condition_column)
    num_columns += 1
    num_rows += 1
    """
    Next we want to add in the columns corresponding to our auxiliary variables
    """
    auxiliary_columns = []
    for i in range(num_rows):
        new_column = [0.0] * num_rows
        new_column[i] = 1.0
        auxiliary_columns.append(new_column)
    columns.extend(auxiliary_columns)
    num_columns += len(auxiliary_columns)
    """
    This then represents A for Ax=b
    """
    A = matrix(columns)
    print columns[2]
    """
    Now we create our objective vector c
    """
    c = matrix([0.0] * (num_columns - len(auxiliary_columns)) + [1.0] * len(auxiliary_columns))
    """
    Next we create b
    """
    b = matrix([0.0] * (num_rows - 1) + [1.0])  # the extra 1 is for the sum condition row
    """
    Next we create G which is just the negative identity
    """
    columns_of_G = []
    for i in range(num_columns):
        new_column = [0.0] * num_columns
        new_column[i] = -1.0
        columns_of_G.append(new_column)
    G = matrix(columns_of_G)
    """
    And finally h
    """
    h = matrix([0.0] * num_columns)
    return c, G, h, A, b


def get_solution(E):
    c, G, h, A, b = retrieve_components(E)
    dims = {
        'l': G.size[0],
        'q': [],
        's': []
    }
    solution = solvers.conelp(c, G, h, dims, A, b)
    return solution


def grab_exact_solution(A, basis):
    # basis should be a list of columns
    # A should be an augmented matrix
    rows_left = list(range(len(A)))
    c_to_r = {}
    for h in range(len(basis)):
        c = basis[h]
        found_a_good_row = False
        for x in range(len(rows_left)):
            r = rows_left[x]
            pivot_el = A[r][c]
            if pivot_el != 0:
                rows_left.pop(x)
                c_to_r[c] = r
                found_a_good_row = True
                break
        if not found_a_good_row:
            continue
        # first we get the column in A to be one
        pivot_el = A[r][c]
        A[r] = [el / pivot_el for el in A[r]]
        # next we zero out this column in all the other rows
        for j in [i for i in range(len(A)) if i != r]:
            multiple = A[j][c]
            A[j] = [A[j][k] - multiple * A[r][k] for k in range(len(A[r]))]
    # now that we've row reduced, we go back through our basis to grab our solution
    solution = [Fraction(0)] * (len(A[0]) - 1)
    for c in c_to_r:
        r = c_to_r[c]
        solution[c] = A[r][-1]
    return solution


def rank(A):
    C = []
    for r in range(len(A)):
        C.append([float(A[r][c]) for c in range(len(A[r]))])
    return matrix_rank(C, 10**-12)


def solve_kirky_with_extension(E):
    RE = reduce_matrix(E)
    sum = 0
    column = []
    for row in RE:
        sum += row[2]
        column.append(row[2])
    print 'SUM: %s' % sum
    print 'COL: %s' % column
    solution = get_solution(RE)
    print solution['primal objective']
    if solution['primal objective'] > 10 ** -5: #TODO better decision making here
        print 'HELLO'
        return
    else:
        # now I want to grab the biggest n where n is the number of rows in E...
        # so we're getting a basis w.r.t to our huge system... So that's what
        # we really need to get at.
        m = len(RE[0])
        print solution['x']
        solution_copy = [x for x in solution['x']][:m]
        n = len(RE)
        maxes = set()
        print 'N: %s' % n
        print 'M: %s' % m
        for i in range(n + 1):  # we have the plus one because
                                # our lp problem has an additional row
            _max = max(solution_copy)
            for j in range(len(solution_copy)):
                if solution_copy[j] == _max:
                    solution_copy.pop(j)
                    break
            #print solution_copy
            maxes.add(_max)
        cleaned_solution = []
        basis = []
        for i in range(m + 1):
            if solution['x'][i] in maxes:
                basis.append(i)
                #cleaned_solution.append(solution['x'][i])
            else:
                pass
                #cleaned_solution.append(0.0)
        # now we add on the extra columns for augmentation and sum condition
        RE = reduce_matrix_exact(E)
        RE = [(RE[i] + [Fraction(0), Fraction(0)]) for i in range(len(RE))]
        print len(basis)
        # and we add the sum condition row
        RE.append([Fraction(1)] * (len(RE[0]) - 2) + [Fraction(-1), Fraction(1)])
        cleaned_solution = grab_exact_solution(RE, basis)[:-1]
        print len(cleaned_solution)
        return cleaned_solution
