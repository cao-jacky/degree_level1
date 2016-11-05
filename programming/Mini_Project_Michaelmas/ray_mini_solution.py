''' mini project solution module 2015 - 2016 '''
import numpy

def refraction_2d (incident_rays, planar_surface):
    ''' example implementation of mini project 2015-2016 
        (NOTE: details of parameters in numpy help format could be added here)'''
    
    # set up array to hold results
    refracted_rays = numpy.zeros(incident_rays.shape, dtype=float)
    # NOTE: this is ONE option to avoid a transpose when constructing the returned array of results
    
    # calculate critical angle
    if planar_surface[5] >= planar_surface[4]:
        # total internal reflection is not possible - save computational time by not calculating the critical angle
        tir_possible = False
    else:
        tir_possible = True
        critical_angle = numpy.arcsin(planar_surface[5]/planar_surface[4])
    
    # calculate angles to surface normal of incident rays
    planar_angle = numpy.arctan((planar_surface[2]-planar_surface[0])/(planar_surface[3]-planar_surface[1]))
    incident_angles =  incident_rays[:,2] - (numpy.pi/2.0) - planar_angle # transform ray angles to be with respect to normal to surface
    
    # handle incident rays exceeding the critical angle
    if tir_possible and(abs(incident_angles) > critical_angle).any():
        raise Exception, "at least one incident ray exceeds the critical angle"
    
    # calculate gradients and intercepts of incoming rays and surface
    ray_gradients = numpy.tan((numpy.pi/2.0) - incident_rays[:,2])
    ray_intercepts = incident_rays[:,1] - (ray_gradients * incident_rays[:, 0])
    surface_gradient = (planar_surface[3]-planar_surface[1])/(planar_surface[2]-planar_surface[0])
    surface_intercept = planar_surface[1] - (surface_gradient * planar_surface[0])

    # calculate points of intersection of rays with surface...
    # horizontal
    refracted_rays[:,0] = (surface_intercept - ray_intercepts) / (ray_gradients - surface_gradient)
    # vertical
    refracted_rays[:,1] = (ray_gradients * refracted_rays[:,0]) + ray_intercepts
    
    # calculate directions of refracted rays
    refracted_angles = numpy.arcsin(planar_surface[4]*numpy.sin(incident_angles)/planar_surface[5]) # Snell's law
    refracted_rays[:,2] = refracted_angles + (numpy.pi/2.0) + planar_angle # transform output ray angles to be clockwise from vertical
    
    return refracted_rays

# test code
if __name__ == "__main__":
    print refraction_2d(numpy.array([[0.0,5.0, numpy.pi/2.0],[0.0,5.0, 1.1*numpy.pi/2.0]]), numpy.array([5.0,2.0,6.0,8.0,1.0,1.33]))