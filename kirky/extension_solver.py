from cvxopt import matrix, solvers

"""
For information on the extension being used see: http://cvxopt.org/userguide/coneprog.html#linear-programming
"""


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


def solve_kirky_with_extension(E):
    solution = get_solution(E)
    if solution['primal objective'] > 10 ** -5: #TODO better decision making here
        return
    else:
        # now I want to grab the biggest n where n is the number of rows in E...
        m = len(E[0])
        print solution['x']
        solution_copy = [x for x in solution['x']][:m]
        n = len(E)
        maxes = set()
        for i in range(n + 1):  # we have the plus one because
                                # our lp problem has an additional row
            _max = max(solution_copy)
            for j in range(len(solution_copy)):
                if solution_copy[j] == _max:
                    solution_copy.pop(j)
                    break
            maxes.add(_max)
        cleaned_solution = []
        for i in range(m):
            if solution['x'][i] in maxes:
                cleaned_solution.append(solution['x'][i])
            else:
                cleaned_solution.append(0.0)
        return cleaned_solution
