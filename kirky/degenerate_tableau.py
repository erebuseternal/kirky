from fractions import Fraction

def add_vectors(v1, v2):
    i = 0
    result = []
    for e1 in v1:
        e2 = v2[i]
        result.append(e1 + e2)
        i += 1
    return result


class DegenerateTableau(object):

    def __init__(self, c, d, A, b):
        # setup the object vector (first row of tableau as written)
        self.objective_vector = [e for e in c]
        # get the objective value (to make the matrix augmented)
        self.objective_value = -d
        # setup the constraint vectors as values
        self.constraint_vectors = []
        self.constraint_values = []
        # we will find the pivot row while we are at it
        for i in range(len(A)):
            constraint_vector = [e for e in A[i][:]]
            self.constraint_vectors.append(constraint_vector)
            self.constraint_values.append(b[i])
        # and finally find the basis
        ###
        self.last_row = None
        ###
        self.basis = []
        self.find_basis()
        self.mend_objective()
        self.pivot_column = None
        self.pivot_row = None
        self.find_pivot()
        self.solution = []
        self.get_solution()

    def find_basis(self, given=None):
        used_i = []
        if given is not None:
            self.basis = [given]
            for i in range(len(self.constraint_vectors)):
                if self.constraint_vectors[i][given] == 1:
                    used_i.append(i)
        else:
            self.basis = []
        for j in range(len(self.objective_vector)):
            if len(self.basis) == len(self.constraint_vectors):
                break
            if j == given:
                continue
            found_one = False
            just_one = True
            for i in range(len(self.constraint_vectors)):
                if i in used_i:
                    continue
                if self.constraint_vectors[i][j] == 1:
                    if not found_one:
                        found_one = True
                        current_i = i
                    else:
                        just_one = False
                elif self.constraint_vectors[i][j] != 0:
                    just_one = False
                    break
            if just_one and found_one:
                self.basis.append(j)
                used_i.append(current_i)

    def mend_objective(self):
        for i in self.basis:
            if self.objective_vector[i] == 0:
                continue
            coefficient = self.objective_vector[i]
            for j in range(0, len(self.constraint_vectors)):
                if self.constraint_vectors[j][i] == 1:
                    break
            for k in range(len(self.constraint_vectors[j])):
                if k == i:
                    continue
                element = self.constraint_vectors[j][k]
                self.objective_vector[k] += coefficient * -1 * element
            value = self.constraint_values[j]
            self.objective_value += -1 * coefficient * value
            # and finally set the non zero element we found that started this whole thing to zero
            self.objective_vector[i] = Fraction(0, 1)

    def find_pivot(self):
        self.pivot_row = None
        self.pivot_column = None
        # because this is the degenerate case we use Bland's rule
        for i in range(len(self.objective_vector)):
            if self.objective_vector[i] < 0:
                self.pivot_column = i
                break
        if self.pivot_column is None:
            return
        smallest_ratio = 10**16
        for i in range(len(self.constraint_vectors)):
            ###
            #if i == self.last_row:
            #    continue
            ###
            column_value = self.constraint_vectors[i][self.pivot_column]
            if column_value <= 0:
                i += 1
                continue
            constraint_value = self.constraint_values[i]
            ratio = constraint_value / column_value
            if ratio < smallest_ratio:
                smallest_ratio = ratio
                self.pivot_row = i

    def get_solution(self):
        self.solution = [Fraction(0, 1)] * len(self.constraint_vectors[0])
        for j in self.basis:
            for i in range(len(self.constraint_vectors)):
                if self.constraint_vectors[i][j] == 1.0:
                    self.solution[j] = self.constraint_values[i]

    def is_solved(self):
        if self.pivot_column is None and self.pivot_row is None:
            return True
        return False

    def is_unbounded(self):
        if self.pivot_column is not None and self.pivot_row is None:
            return True
        return False

    def pivot(self):
        divisor = self.constraint_vectors[self.pivot_row][self.pivot_column]
        pivot_vector = [e / divisor for e in self.constraint_vectors[self.pivot_row]]
        self.constraint_vectors[self.pivot_row] = pivot_vector
        pivot_value = self.constraint_values[self.pivot_row] / divisor
        self.constraint_values[self.pivot_row] = pivot_value
        for i in range(len(self.constraint_vectors)):
            vector = self.constraint_vectors[i]
            if vector == pivot_vector:
                continue
            coefficient = -1 * vector[self.pivot_column]
            adjustment_vector = [e * coefficient for e in pivot_vector]
            new_vector = add_vectors(vector, adjustment_vector)
            self.constraint_vectors[i] = new_vector
            self.constraint_values[i] += pivot_value * coefficient
        coefficient = -1 * self.objective_vector[self.pivot_column]
        adjustment_vector = [e * coefficient for e in pivot_vector]
        self.objective_vector = add_vectors(self.objective_vector, adjustment_vector)
        self.objective_value += pivot_value * coefficient
        # last but not least we find the new basis and pivot
        self.find_basis(self.pivot_column)
        ###
        self.last_row = self.pivot_row
        ###
        self.find_pivot()
        self.get_solution()


def solve_kirky(E):
    num_weights = len(E[0])
    # we begin by creating the row that will represent the positive sum condition
    sum_condition_row = [Fraction(1)] * num_weights + [Fraction(-1)]
    sum_condition_value = Fraction(1)
    # we now add on the extra column for all the rows of E
    for row in E:
        row.append(Fraction(0))
    # we now augment E with the sum condition row
    E.append(sum_condition_row)
    b = [Fraction(0)] * (len(E) - 1) + [sum_condition_value]
    # now we will use a phase I approach for simplex to find a feasible solution
    # if it exists
    # we begin by adding in our auxiliary variables
    num_auxiliary_variables = len(E)
    i = 0
    for row in E:
        addition = [Fraction(0)] * num_auxiliary_variables
        addition[i] = Fraction(1)
        row.extend(addition)
        i += 1
    # we create our objective function
    c = [Fraction(0)] * (len(E[0]) - num_auxiliary_variables) + [Fraction(1)] * num_auxiliary_variables
    d = Fraction(0)
    # now we solve this phase I step
    tabby = DegenerateTableau(c, d, E, b)
    while not tabby.is_solved():
        tabby.pivot()
    if tabby.objective_value != 0:
        return None
    else:
        tabby.get_solution()
        return tabby.solution[:num_weights]




