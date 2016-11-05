import numpy
import matplotlib.pyplot as pyplot

def plotfn(x, xlabel, ylabel, title):
    #The function which is being used
    x1 = x
    y = (2 * (x1 ** 2)) + (3 * x1) + 4

    #Plotting function
    pyplot.plot(x,y)
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    pyplot.title(title)

    pyplot.show()

    return plotfn

if __name__ == '__main__':

    import ex12_2
    import numpy

    ex12_2.plotfn(numpy.linspace(0,9,10),'test x', 'testy', 'lel')
