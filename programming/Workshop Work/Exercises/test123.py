'''
Program to tabulate some integers, their squares and their cubes
'''
import numpy

size = 6
numbers = numpy.zeros(size, dtype = float)
squares = numpy.zeros(size, dtype = float)
cubes = numpy.zeros(size, dtype = float)
for i in range(size):
    numbers[i] = i
    squares[i] = i**2
    cubes[i] = i**3

print numbers
print squares
print cubes
