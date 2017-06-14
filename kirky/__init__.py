import numpy as np
from .block import Edge, Frame
from .draw import DrawEdge
from pyx import canvas
from .helpers import common_denominator
from fractions import Fraction, gcd
from .tableau import solve_kirky
from .extension_solver import solve_kirky_with_extension


class Kirchhoff(object):

    def __init__(self, matrix):
        """
        Input:
            matrix - of the matrix we are wanting to generate a Kirchhoff graph for
                    we should have put it in [IB] form and the matrix input should be
                    the B part with elements that are python Fractions

        During initiation of the class we generate the frame that will be used during
        the course of the algorithm.

        NOTES:
            (a) we determine the number of dimensions our frame will have
            (b) the number of vectors we will have is equal to the number of coordinate vectors
                (given by the dimensions of the frame) plus the number of cross vectors
                (given by the number of columns in our matrix)
            (c) we determine how closely the vertices in the frame will be to each other
                along each dimension. This gets saved in steps as the number of vertices
                in [0, 1) along that dimension. So if we have a distance of 1/2
                between vertices along dimension 0 then step[0] = 2.
            (d) we get the shape of the smallest possible frame that contains all the vectors
                we have
            (e) we create a frame of the appropriate shape filled only with coordinate vectors
            (f) the next three lines fill the frame with all the possible cross vectors
        """
        self.matrix = matrix
        self.dimensions = matrix.shape[0]                                                               # (a)
        self.num_vectors = self.dimensions + matrix.shape[1]                                            # (b)
        self.steps = self.get_steps()                                                                   # (c)
        self.frame = Frame(self.dimensions, self.num_vectors, self.steps)
        first_frame_shape = self.get_first_frame_shape()                                                # (d)
        self.frame.seed_frame(first_frame_shape)                                                        # (e)
        cross_vectors = self.get_cross_vectors()                                                        # (f)
        for vector in cross_vectors:
            self.frame.populate(vector)

    def get_steps(self):
        """
        For each dimension we get the number of steps to take in moving along vertices from 0
        to 1 along that dimension. If along that dimension our vectors have components 1/2 1
        2/3 0, then we can see that the smallest distance we can move in that dimension by using
        these vectors is 1/6. Therefore we will have 6 steps between 0 and 1 to get all the possible
        vertices along this dimension. We essentially do this for each dimension and output the
        results in a list.
        """
        steps = []
        for row in self.matrix:
            steps.append(common_denominator([fraction for fraction in row]))
        return steps

    def get_first_frame_shape(self):
        """
        This returns the smallest shape that still gives us at least one of each vector.
        """
        shape = []
        for row in self.matrix:
            max_element = max([abs(element) for element in row] + [Fraction(1)])
            shape.append(max_element)
        return shape

    def get_cross_vectors(self):
        """
        This creates an edge object for each cross vector and returns them in a list.

        NOTES:
            (a) each column in our matrix is a cross vector
            (b) the id's are basically the vector label
        """
        edges = []
        trans = np.transpose(self.matrix)                                                               # (a)
        id = self.dimensions                                                                            # (b)
        for row in trans:
            head = [element for element in row]
            tail = [Fraction(0) for _ in range(self.dimensions)]
            edge = Edge(tail, head, id)
            edges.append(edge)
            id += 1
        return edges

    def generate_linear_system(self):
        """
        Refer to the report, but this generates the rows of our constraint matrix E which
        if solved will give us the weightings of the frame vectors which make the result
        Kirchhoff
        """
        rows = []
        for position in self.frame.vertices:
            vertex = self.frame.vertices[position]
            for i in range(self.matrix.shape[1]):
                row = [Fraction(0) for _ in range(len(self.frame.coordinate_vectors) + len(self.frame.cross_vectors))]
                for j in range(self.dimensions):
                    multiplier = self.matrix[j, i]
                    in_edge, out_edge = vertex.cut[j]
                    if in_edge:
                        row[in_edge.pin] = -1 * multiplier
                    if out_edge:
                        row[out_edge.pin] = 1 * multiplier
                in_edge, out_edge = vertex.cut[self.dimensions + i]
                if in_edge:
                    row[in_edge.pin] = Fraction(1)
                if out_edge:
                    row[out_edge.pin] = Fraction(-1)
                rows.append(row)
        return rows

    def solve(self, linear_system):
        """
        takes a linear system and returns a solution if one exists (the solution will be None
        if a solution could not be found)
        """
        vector_solution = solve_kirky(linear_system)
        return vector_solution

    def normalize_solution(self, vector_solution):
        """
        This takes the solution (composed of fractions) and multiplies through by a single constant
        (which is found during the process of this method) to generate a fully integer solution.
        """
        denominators = {weight.denominator for weight in vector_solution}
        while len(denominators) > 1:
            l = list(denominators)
            first, second = l[:2]
            new_denominator = (first / gcd(first, second)) * second
            denominators.remove(first)
            denominators.remove(second)
            denominators.add(new_denominator)
        common_denominator = list(denominators)[0]
        normalized_solution = [weight * common_denominator for weight in vector_solution]
        return normalized_solution

    def set_edge_weights(self, solution):
        """
        This takes a solution and sets the corresponding edge weights in the frame. The importnat thing
        to note here is that the solution is sorted by the pin numbers of the edges.
        """
        edges = [edge for edge in self.frame.coordinate_vectors] + [edge for edge in self.frame.cross_vectors]
        edges = sorted(edges, key=lambda edge: edge.pin)
        print len(edges)
        print len(solution)
        for i in range(len(edges)):
            edges[i].weight = solution[i]

    def draw_edges(self, canvas):
        count = 1
        print 'DRAWING COORDINATE VECTORS'
        for edge in self.frame.coordinate_vectors:
            print 'drawing edge %s/%s' % (count, len(self.frame.coordinate_vectors))
            DrawEdge(edge, canvas)
            count += 1
        count = 1
        print 'DRAWING CROSS VECTORS'
        for edge in self.frame.cross_vectors:
            print 'drawing edge %s/%s' % (count, len(self.frame.cross_vectors))
            DrawEdge(edge, canvas)
            count += 1

    def find(self, file='', exact=True):
        """
        This is the meat of the algorithm. It takes the frame and runs the steps outlined
        in the report. So if you've read the report this should be pretty self explanatory.
        """
        dimension = 0
        while True:
            print '--> generating linear system'
            linear_system = self.generate_linear_system()
            print '     size is %s, %s' % (len(linear_system), len(linear_system[0]))
            print '--> trying to find aa solution'
            if exact:
            	vector_solution = self.solve(linear_system)
            else:
                vector_solution = solve_kirky_with_extension(linear_system)
            if vector_solution is None:
                print '--> solution not found'
                print '--> doubling along dimension %s' % dimension
                self.frame.double(dimension)
                dimension = (dimension + 1) % self.dimensions
            else:
                break
        if exact:
            vector_solution = self.normalize_solution(vector_solution)
        self.set_edge_weights(vector_solution)
        if file:
            c = canvas.canvas()
            self.draw_edges(c)
            c.writePDFfile(file)

    def incidence(self):
        edgetionary = {}
        incidence_matrix = []
        positions = []
        label_matrix = []
        edges = [edge for edge in self.frame.coordinate_vectors] + [edge for edge in self.frame.coordinate_vectors]
        counter = 0
        for edge in edges:
            if edge.weight != 0:
                edgetionary[edge.pin] = counter
                row = [Fraction(0)] * self.num_vectors
                row[edge.id] = Fraction(1)
                label_matrix.append(row)
                counter += 1
        num_columns = len(edgetionary)
        for position in self.frame.vertices:
            non_zero_cut = False
            row = [Fraction(0)] * num_columns
            for double in self.frame.vertices[position].cut:
                if double[0] is not None and double[0].pin in edgetionary:
                    column = edgetionary[double[0].pin]
                    row[column] = -double[0].weight
                    non_zero_cut = True
                if double[1] is not None and double[1].pin in edgetionary:
                    column = edgetionary[double[1].pin]
                    row[column] = double[1].weight
                    non_zero_cut = True
            if non_zero_cut:
                incidence_matrix.append(row)
                positions.append(list(position))
        return incidence_matrix, label_matrix, positions

    def graphical_indicence(self):
        incidence_matrix = self.incidence()[0]
        graphical = []
        for row in incidence_matrix:
            graphical_row = []
            for e in row:
                if e == 0:
                    graphical_row.append(0)
                else:
                    graphical_row.append(int(e / abs(e)))
            graphical.append(graphical_row)
        return graphical

    def see_frame(self, file):
        for edge in self.frame.coordinate_vectors:
            edge.weight = edge.pin
        for edge in self.frame.cross_vectors:
            edge.weight = edge.pin
        c = canvas.canvas()
        self.draw_edges(c)
        c.writePDFfile(file)

    def try_size(self, frame_shape, file=''):
        # first we setup the block anew
        self.frame = Frame(self.dimensions, self.num_vectors, self.steps)
        first_frame_shape = self.get_first_frame_shape()
        self.frame.seed_frame(first_frame_shape)
        cross_vectors = self.get_cross_vectors()
        for vector in cross_vectors:
            self.frame.populate(vector)
        # now we try to grow our block to the right size
        if not self.frame.grow_to_size(frame_shape):
            return None
        print '--> generating linear system'
        linear_system = self.generate_linear_system()
        print '     size is %s, %s' % (len(linear_system), len(linear_system[0]))
        print '--> trying to find a solution'
        vector_solution = self.solve(linear_system)
        if vector_solution is None:
            print '--> solution not found'
            return None
        else:
            print '--> solution found'
        solution = self.normalize_solution(vector_solution)
        self.set_edge_weights(solution)
        if file:
            c = canvas.canvas()
            self.draw_edges(c)
            c.writePDFfile(file)





