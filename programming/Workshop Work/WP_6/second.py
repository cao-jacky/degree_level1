import numpy
import matplotlib.pyplot as pyplot

#Function
x = numpy.linspace(0,9.0,10)
y = 2 * x

pyplot.plot(x,y)

#Axis labels
pyplot.xlabel('x value in radians')
pyplot.ylabel('f(x)')
pyplot.title('y=f(x)')

pyplot.show()
