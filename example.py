from kirky import Kirchhoff
import numpy as np
from fractions import Fraction

"""
To generate the Kirchhoff Graph for the matrix
1 0 2 1
0 1 1 2
"""
matrix = np.array([[Fraction(2),Fraction(1)], [Fraction(1),Fraction(2)]])
k = Kirchhoff(matrix)
k.find('drawings/kirchhoff1')

"""
To generate the Kirchhoff Graph for the matrix
1 0 -1 1 1
0 1  1 2 1
"""
#matrix = np.array([[Fraction(-1),Fraction(1),Fraction(1)],[Fraction(1),Fraction(2),Fraction(1)]])
#k = Kirchhoff(matrix)
#k.find('drawings/kirchhoff2')

matrix = np.array([[Fraction(1)], [Fraction(1)]])
k = Kirchhoff(matrix)
k.find('drawings/kirchhoff2')

"""
To generate the Kirchhoff Graph for the matrix
1 0 -1 1  3
0 1  1 2 -1
"""
matrix = np.array([[Fraction(-1),Fraction(1), Fraction(3)],[Fraction(1),Fraction(2), Fraction(-1)]])
k = Kirchhoff(matrix)
k.find('drawings/kirchhoff3')

"""
To generate the Kirchhoff Graph for the matrix
1 0 -3/2 1/2
0 1  1/2 1/2
"""
#matrix = np.array([[Fraction(-3,2),Fraction(1,2)], [Fraction(1,2),Fraction(1,2)]])
#k = Kirchhoff(matrix)
#k.find('drawings/kirchhoff4')
