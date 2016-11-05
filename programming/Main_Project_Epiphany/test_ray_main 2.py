''' simple test harness for ray_main.py main project 2015 - 2016. Simulates marker program importing and calling
submitted mini project. This DOES NOT test the correctness of the returned results. It just checks that the ray_main module can be imported, that the various task functions can be called with example arguments, and it prints the contents and shape of the returned results'''
import numpy
import ray_main

print " **** IMPORTANT: this test program MUST be run from a command window with:"
print " **** python test_ray_main.py"
print " **** it must be run from the same folder as ray_main.py"
print "---------------------------"

# refraction_2d (or equivalent) is still required in the main project, but will not be marked directly.
print "refraction_2d (from the mini project)"
rr = ray_main.refraction_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]), numpy.array([4.5, 2.0, 5.5, 8.0,1.0,1.33]))
print rr
print "returned data type is ",type(rr), " -- this should be  <type 'numpy.ndarray'>"
print "returned data shape is ",rr.shape, " -- this should be (2, 3) or (2L, 3L) for this example of 2 incident rays"
print "---------------------------"

print "refraction_2d_sph (TASK 1)"
rr = ray_main.refraction_2d_sph(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]), numpy.array([0.0, 5.0,0.0,-5.0,1.0,1.33,15.0]))
print rr
print "returned data type is ",type(rr), " -- this should be  <type 'numpy.ndarray'>"
print "returned data shape is ",rr.shape, " -- this should be (2, 3) or (2L, 3L) for this example of 2 incident rays"
print "---------------------------"

print "refraction_2d_DET (TASK 2)"
rr = ray_main.refraction_2d_det(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]), 5.0)
print rr
print "returned data type is ",type(rr), " -- this should be  <type 'numpy.ndarray'>"
print "returned data shape is ",rr.shape, " -- this should be (2, 3) or (2L, 3L) for this example of 2 incident rays"
print "---------------------------"

print "trace_2d (TASK 3)"
rrp = ray_main.trace_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]),
                            [["PLA", numpy.array([4.5, 2.0, 5.5, 8.0,1.0,1.33])],
                             ["DET", 10.0]])
print rrp
print "returned data type is ",type(rrp), " -- this should be  <type 'numpy.ndarray'>"
print "returned data shape is ",rrp.shape, " -- this should be (2, 2, 3) or (2L, 2L, 3L) for this example of 2 surfaces and 2 incident rays"
print "---------------------------"

print "plot_trace_2d (TASK 4) (this test depends on the TASK 3 test functioning correctly)"
ray_main.plot_trace_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]),rrp,
                        [["PLA", numpy.array([4.5, 2.0, 5.5, 8.0,1.0,1.33])],
                         ["DET", 10.0]])
print "---------------------------"

print "evaluate_trace_2d (TASK 5) (this test depends on the TASK 3 test functioning correctly)"
frac = ray_main.evaluate_trace_2d(rrp, 0.6941) # 0.1 is not guaranteed to be a sensible value - adjust as required to test
print frac
print "returned data type is ",type(frac), " -- this should be", type(0.6941)
print "---------------------------"

print "optimise_surf_rad_2d (TASK 6)"
rad_opt = ray_main.optimize_surf_rad_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]),
                                 [["PLA", numpy.array([4.5, -8.0, 5.5, 8.0,1.0,1.33])],
                                  ["SPH", numpy.array([0.0, -10.0,0.0,10.0,1.0,1.33,15.0])],
                                  ["DET", 15.0]], 0.1, numpy.array([1]))
print rad_opt
print "returned data type is ",type(rad_opt), " -- this should be", type(numpy.array([1.0]))
print "returned data shape is ",rad_opt.shape, " -- this should be", numpy.array([1.0]).shape

print "---------------------------"
