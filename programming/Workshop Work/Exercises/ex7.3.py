"""lets create some code to generate the first 10 elements of Fibonacci series of integers"""

import numpy
x = numpy.array([0,1])
y = numpy.array([1])

next_value = 0

for i in range(len(x)):
    next_value = next_value + x[i]

print x

#more testing
size = 10
billyray = numpy.zeros(size, dtype=float) #billyray is an index
billyray[0] = 0

for i in range(size):
    billyray[i] = billyray[i] + billyray[i]

print billyray[i]
