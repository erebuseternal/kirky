from .tableau import Simplex
from fractions import Fraction
from time import sleep
from random import sample


def dot(vector1, vector2):
    return sum([vector1[i] * vector2[i] for i in range(len(vector1))])


def negatives(vectors):
    vectors = {tuple(vector) for vector in vectors}
    for vector in vectors:
        reflection = tuple([-e for e in vector])
        if reflection in vectors:
            return True
    return False


class Splitter(object):

    def __init__(self, vectors, exact=False):
        vectors = self.clean(vectors)
        self.infeasible = False
        if len(vectors) == 0:
            print 'no vectors found'
            self.infeasible = True
        else:
            self.exact = exact
            self.A = []
            self.dimension = len(vectors[0])
            for i in range(len(vectors)):
                if not self.exact:
                    row = [e for e in vectors[i]] + [-e for e in vectors[i]] + [-1.0]
                    id_row = [0.0 for j in range(len(vectors))]
                    id_row[i] = -1.0
                    row += id_row
                else:
                    row = [e for e in vectors[i]] + [-e for e in vectors[i]] + [Fraction(-1, 1)]
                    id_row = [Fraction(0, 1) for j in range(len(vectors))]
                    id_row[i] = Fraction(-1, 1)
                    row += id_row
                self.A.append(row)
            if self.exact:
                self.b = [Fraction(1, 1)] * len(self.A)
            else:
                self.b = [1.0] * len(self.A)
            self.simplex = Simplex(self.A, self.b, [0.0] * len(self.A[0]), exact=self.exact)

    def clean(self, vectors):
        good_vectors = []
        for vector in vectors:
            if dot(vector, vector) == 0:
                continue
            good_vectors.append(vector)
        return good_vectors

    def get_simplex_solution(self):
        self.simplex.prep()
        if not self.simplex.is_feasible:
            return None
        else:
            return self.simplex.solution

    def convert_simplex_solution(self, simplex_solution):
        norm = []
        for i in range(self.dimension):
            w = simplex_solution[i] - simplex_solution[i + self.dimension]
            norm.append(w)
        return norm

    def split(self):
        if self.infeasible:
            return None
        simplex_solution = self.get_simplex_solution()
        if simplex_solution is None:
            return None
        else:
            return self.convert_simplex_solution(simplex_solution)

