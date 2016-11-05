import numpy
stuffs = numpy.arange(1,101)
stuffs.shape = (10,10)
print stuffs

#we are now going to slice the 5's column
slice = stuffs[0:10, 4]
print slice

#slicing an array from 35-38, 45-48, 55-58
slice = stuffs [3:6, 4:8]
print slice

"""
#flattening the above matrix
woopty = slice
#a.flatten(woopty)
numpy.ndarray.flatten(woopty)
print woopty

#IT DOESN'T FUCKING WORK AND THER'S NO BLOODY INTERNET SO WE CANNOT FUCKING CHECK ALSO FACEBOOK DOESN'T WORK YOU CUNTS.
numpy.array(woopty)
numpy.ndarray.flatten(woopty)
woopty.flatten()
print woopty
"""
import numpy as np

fuck = slice
a = np.array(fuck)
a = a.flatten()
print a

"""
#test again
a = np.reshape(a (1,np.product(a.shape)))
print a
"""
