import unittest
from numpy.linalg import matrix_rank
from numpy import array
from extension_solver import *

class TestSolverSteps(unittest.TestCase):

    def test_reduce_matrix(self):
        A1 = [
	    [1, 0, 1, 0],
            [1, 0, 1, 0],
            [0, 0, 0, 1]  
	]
        A2 = [
            [1, 0, 1],
            [1, 1, 0],
            [0, 0, 1],
            [0, 0, 1]
        ]
        A3 = [
            [0,     100],
            [0,     237],
            [10**-3, 0] 
        ]
        A1 = reduce_matrix(A1)
	self.assertEqual(matrix_rank(array(A1)), 2)
	self.assertEqual(len(A1), 2)
	self.assertEqual(len(A1[0]), 4)
	A2 = reduce_matrix(A2)
        self.assertEqual(matrix_rank(array(A2)), 3)
        self.assertEqual(len(A2), 3)
	self.assertEqual(len(A2), 3)  
	A3 = reduce_matrix(A3)
        self.assertEqual(matrix_rank(array(A3)), 2)
        self.assertEqual(len(A3), 2)
	self.assertEqual(len(A3), 2) 
    
    def test_create_augmented_matrix(self):
        E = [
            [1, 2, 3],
            [4, 5, 6]
        ]
        A = create_augmented_matrix(E)
        expected_A = [
            [1., 2., 3., 0.],
            [4., 5., 6., 0.],
            [1., 1., 1., -1.],
        ]
	self.assertEqual(A, expected_A)
