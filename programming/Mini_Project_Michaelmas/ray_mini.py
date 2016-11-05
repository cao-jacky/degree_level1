"""A module dedicated to evaluating horizontal and vertical positions of the
points of intersection between an array of incident rays and a specified
surface. Also evaluates the clockwise angles (in radians) with repsect to the
vertical of the refracted rays."""

USER = "Jacky Cao"
USER_ID = "bbvw84"

import numpy

def refraction_2d (incident_rays, planar_surface):
    '''A function which can calculate the location of a refracted ray and also
    the angle the refracted ray makes with the vertical. Input requires is an
    array for both incident rays and planar surface.'''

    #Working out the angle the plane is tilted at:
    vertical_1 = planar_surface[0] - planar_surface[2]
    horizontal_1 = planar_surface[1] - planar_surface[3]
    angle_plane = numpy.arctan((vertical_1/horizontal_1))
    #print ap

    #gradients
    #plane
    gradient_plane = (planar_surface[3] - planar_surface[1])/(planar_surface[2] - planar_surface[0])
    #incident ray
    gradient_ray = numpy.tan((numpy.pi / 2) - incident_rays[:,2])

    #y-intercepts
    #planar equation
    intercept_plane = (- planar_surface[0] * gradient_plane) + planar_surface[1]
    #incident ray
    intercept_ray = (- incident_rays[:,0] * gradient_ray) + incident_rays[:,1]

    #point of intersection with plane
    #x-coord
    x_intersect_plane = (intercept_plane - intercept_ray)/(gradient_ray - gradient_plane)
    #y-coord
    y_intersect_plane = ((intercept_ray * gradient_plane) - (intercept_plane * gradient_ray)) / (gradient_plane - gradient_ray)

    #angle of incidence - corrected
    angle_incidence = incident_rays[:,2] - (numpy.pi / 2) - angle_plane

    numpy.seterr(all = 'ignore') #ignores those pesky runtime error
    critical_angle = numpy.arcsin((planar_surface[5] / planar_surface[4]))

    #Critical angle function - corrected
    if planar_surface[5] >= planar_surface[4]:
        tir_possible = False
    else:
        tir_possible = True

    #function terminator if the angles given are greater than the critical angle
    if tir_possible and any(angle_incidence) >= critical_angle:
        raise Exception, "at least one incident ray exceeds the critical angle"
        return
    else:
        print "Incident angles are valid, function will proceed."

    #Corrected code for the angle of refraction
    refracted_angles = numpy.arcsin(planar_surface[4]*numpy.sin(angle_incidence)/planar_surface[5]) # Snell's law
    refracted_angles = refracted_angles + (numpy.pi/2.0) + angle_plane # transform output ray angles to be clockwise from vertical

    #Outputting final array for refracted rays
    final_ray = numpy.array([x_intersect_plane, y_intersect_plane, refracted_angles]) #the rays are shown vertically
    final_ray = numpy.rot90(final_ray, 3) #next two lines of code are to make it look the same as the initial incident rays array
    refracted_rays = numpy.fliplr(final_ray)

    """
    print "refracted_rays:"
    print refracted_rays
    """

    return refracted_rays

#Testing code here
if __name__ =='__main__':

    print "beginning of test" #Useful to show me where code is actually being tested and if anything from the function has 'leaked' out

    import ray_mini #should be a self contained module, apart from the inital data - as shown below
    import numpy

    #test data
    #incident_rays = numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]])
    #planar_surface = numpy.array([5.0,2.0,6.0,8.0,1.0,1.33])

    ray_mini.refraction_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]), numpy.array([5.0, 2.0, 6.0, 8.0,1.0,1.33]))

    #incident_rays = numpy.array([[1.0,0.0,numpy.pi],[1.0,0.0,numpy.pi]]) #more test data
    #planar_surface = numpy.array([0.0,0.0,3.0,0.0,1.33,1.0])

    #print "Printing incident_rays"
    #print incident_rays

    #ray_mini.refraction_2d(incident_rays, planar_surface) #calling the function from the imported module

    print "end of test"
