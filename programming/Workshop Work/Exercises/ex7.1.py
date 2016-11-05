"""
#let's try and create an initial list of integers from 0.0 to 5.0
test = range(0, 5, 1)
print test

#test code because why not
for number in range(10):
    print "The number:", number, "It's square:", number**2
"""

#let's try and do the actual exercise based off the notes
import numpy
# x = numpy.array([range(0,5,1)]) - old code

x = numpy.array([0.,1.,2.,3.,4.,5.]) #creating an array from a list
y = numpy.zeros(len(x), dtype=float) #creating an array to hold the squares
z = numpy.zeros(len(x), dtype=float) #stored cubed values

for i in range(len(x)): #squared x
    y[i] = x[i]**2

for j in range(len(x)): #cubed x
    z[j] = x[j]**3

print x
print y
print z

"""alternatively:

#let's try and do the actual exercise based off the notes
import numpy
# x = numpy.array([range(0,5,1)]) - old code

x = numpy.array([0,1,2,3,4,5]) #creating an array from a list
y = numpy.zeros(len(x), dtype=int) #creating an array to hold the squares
z = numpy.zeros(len(x), dtype=int) #stored cubed values

for i in range(len(x)): #squared x
    y[i] = x[i]**2

for j in range(len(x)): #cubed x
    z[j] = x[j]**3

print x
print y
print z

"""
