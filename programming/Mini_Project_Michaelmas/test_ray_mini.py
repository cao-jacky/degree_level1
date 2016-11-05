''' simple test harness for ray_mini.py mini project 2015 - 2016. Simulates marker program importing and calling
submitted mini project. This DOES NOT test the correctness of the returned results. It just checks that the ray_mini module can be imported, that the refraction_2d function can be called with example arguments, and it prints the contents and shape of the returned array'''
import numpy
import ray_mini

rr = ray_mini.refraction_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]), numpy.array([4.5, 2.0, 5.5, 8.0,1.0,1.33]))
print rr
print "returned data type is ",type(rr), " -- this should be  <type 'numpy.ndarray'>"
print "returned data shape is ",rr.shape, " -- this should be (2, 3) for this example of 2 incident rays"
