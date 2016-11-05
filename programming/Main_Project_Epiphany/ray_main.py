"""An optical ray tracing module with a function for dealing with spherical
surfaces as well as planar surfaces, and also with opaque detector surfaces.
The module is also capable of plotting ray paths and surfaces, together with an
illustration of the distribution of rays at the output. Function can be used to
evaluate the image quality of a give optical design, and a function for
optimising the image quality by automatically adjusting radii of curvature until
maximum image quality is found."""

USER = "Jacky Cao"
USER_ID = "bbvw84"

import numpy
import matplotlib.pyplot as pyplot
import matplotlib.path as pypath
import matplotlib.patches as patches
import matplotlib.hatch as hatch
import sys
import itertools

from scipy.optimize import minimize

#Mini-project original code
def refraction_2d (incident_rays, planar_surface):
    '''Function to calculate the location of a ray (or multiple rays)
    after it has been refracted through a planar surface.

    Parameters:
        incident_rays: a 2-dimensional numpy array of floats: first index is the
        individual incident light ray, second index for each ray is, the horizontal
        coordinate (metres) of the point of origin of the ray, the corresponding
        vertical coordinate, and the clockwise angle (radians) with respect to the
        positive vertical direction of the direction of the incident ray.

        planar_surface: a 1-dimensional numpy array of floats: the horizontal
        coordinate (metres) of one end of the line, the corresponding vertical
        cordinate, the horizontal coordinate of the other end of the line, the
        corresponding vertical coordinate, the refraction index of the medium on
        the incident side of the surface, and the refractive index of the medium
        of the surface.

    Output:
        refracted_rays: a 2-dimensional numpy array of floats: individual refracted
        light ray specified by first index, following values specified by second
        index for each ray; horizontal coordinate (metres) of the point of
        intersection with the planar surface, corresponding vertical coordinate,
        and the clockwise angle (radians) with respect to the positive vertical
        direction of the direction of the refracted ray.'''

    #Working out the angle the plane is tilted at:
    angle_plane = numpy.arctan(((planar_surface[0] - planar_surface[2])/(planar_surface[1] - planar_surface[3])))
    #print angle_plane

    #calcualting the gradients of the various "lines"
    #planar surface
    gradient_plane = (planar_surface[3] - planar_surface[1])/(planar_surface[2] - planar_surface[0])
    #incident ray
    gradient_ray = numpy.tan((numpy.pi / 2) - incident_rays[:,2])

    #calculating the y-intercepts of the various "lines"
    #planar surface
    intercept_plane = (- planar_surface[0] * gradient_plane) + planar_surface[1]
    #incident ray
    intercept_ray = (- incident_rays[:,0] * gradient_ray) + incident_rays[:,1]

    #point of intersection of the ray with plane
    #x-coord
    x_intersect_plane = (intercept_plane - intercept_ray)/(gradient_ray - gradient_plane)
    #y-coord
    y_intersect_plane = ((intercept_ray * gradient_plane) - (intercept_plane * gradient_ray)) / (gradient_plane - gradient_ray)

    #angle of incidence
    angle_incidence = incident_rays[:,2] - (numpy.pi / 2) - angle_plane

    numpy.seterr(all = 'ignore') #ignores those pesky runtime error
    critical_angle = numpy.arcsin((planar_surface[5] / planar_surface[4]))

    #Critical angle function
    if planar_surface[5] >= planar_surface[4]:
        tir_possible = False
    else:
        tir_possible = True

    #function terminator if the angles given are greater than the critical angle
    if tir_possible and any(angle_incidence) >= critical_angle:
        raise Exception, "at least one incident ray exceeds the critical angle"
        return
    #else:
        #print "Incident angles are valid, function will proceed."

    #Corrected code for the angle of refraction
    refracted_angles = numpy.arcsin(planar_surface[4]*numpy.sin(angle_incidence)/planar_surface[5]) # Snell's law
    refracted_angles = refracted_angles + (numpy.pi/2.0) + angle_plane # transform output ray angles to be clockwise from vertical

    #Outputting final array for refracted rays
    final_ray = numpy.array([x_intersect_plane, y_intersect_plane, refracted_angles]) #the rays are shown vertically
    final_ray = numpy.rot90(final_ray, 3) #next two lines of code are to make it look the same as the initial incident rays array
    refracted_rays = numpy.fliplr(final_ray)


    """
    print "refracted_rays: arse"
    print refracted_rays
    """

    return refracted_rays

#Task 1
def refraction_2d_sph (incident_rays, spherical_surface):
    '''Function to calculate the location of a ray (or multiple rays) after it
    has been refracted through a spherical surface.

    Parameters:
        incident_rays: a 2-dimensional numpy array of floats: first index is the
        individual incident light ray, second index for each ray is, the horizontal
        coordinate (metres) of the point of origin of the ray, the corresponding
        vertical coordinate, and the clockwise angle (radians) with respect to the
        positive vertical direction of the direction of the incident ray.

        spherical_surface: a 1-dimensional numpy array of floats: the horizontal
        coordinate (metres) of one end of the line, the corresponding vertical
        cordinate, the horizontal coordinate of the other end of the line, the
        corresponding vertical coordinate, the refraction index of the medium on t
        he incident side of the surface, the refractive index of the medium of the
        surface, the radius of curvature of the surface (metres) - positive radius
        of curvature means a convex surface, negative radius of curvature means a
        concave surface.

    Output:
        refracted_rays: a 2-dimensional numpy array of floats: individual refracted
        light ray specified by first index, following values specified by second
        index for each ray; horizontal coordinate (metres) of the point of
        intersection with the planar surface, corresponding vertical coordinate,
        and the clockwise angle (radians) with respect to the positive vertical
        direction of the direction of the refracted ray.
        '''

    #calculating the gradient of the incident ray
    gradient_ray = numpy.tan((numpy.pi / 2.0) - incident_rays[:,2])
    #print "gradient_ray", gradient_ray

    #calculating the y-intercept of the incident ray
    intercept_ray = (incident_rays[:,1] - (gradient_ray * incident_rays[:, 0]))
    #print "intercept_ray", intercept_ray

    #numpy.seterr(all = 'ignore') #ignores those pesky runtime error

    #Calculating the centre of the circle
    if spherical_surface[6] > 0:
        circle_centre_c = spherical_surface[2] + numpy.sqrt((
                        spherical_surface[6] ** 2) - (
                        spherical_surface[1] ** 2))
        refraction_2d_sph.ccc = circle_centre_c
        ccc = circle_centre_c
        #print "circle_centre_c", circle_centre_c
    else:
        ss = spherical_surface #cutting down on code length
        circle_centre_c = ss[2] + (ss[6] / numpy.abs(ss[6])) * numpy.sqrt(((ss[6] ** 2) - (ss[1] ** 2)))
        refraction_2d_sph.ccc = circle_centre_c
        ccc = circle_centre_c
        #print "circle_centre_cunt", circle_centre_c

    #Working out the point of intersection of the ray with the spherical surface
    #Calculating the x coordinate with a quadratic equation
    ##Various components for quadratic formula
    a = (1 + (gradient_ray ** 2))
    #print "a", a
    b = (2 * intercept_ray * gradient_ray) - (2 * ccc)
    #print "b", b
    c = (ccc ** 2) + (intercept_ray ** 2) - (spherical_surface[6] ** 2)
    #print "c", c

    #Performing the quadratic formula - using an if statement for a concave or convex lens
    if spherical_surface[6] > 0:
        intercept_plane_x = (- b - numpy.sqrt((b ** 2) - (4 * a * c))) / (2 * a)
        #print "intercept_plane_x", intercept_plane_x
    else:
        intercept_plane_x = (- b + numpy.sqrt((b ** 2) + (4 * a * c))) / (2 * a)

    #y coordiante
    intercept_plane_y = (intercept_plane_x * gradient_ray) + intercept_ray
    #print "intercept_plane_y", intercept_plane_y

    #Angle of refraction calculations
    ##Working out a surface normal to the circle at point of intersection

    ###Gradient from point of intersection to centre of the circle
    gradient_radius =  (0 - intercept_plane_y) / (ccc - intercept_plane_x)
    #print "gradient_radius", gradient_radius

    ###Gradient normal
    gradient_normal = - 1 / gradient_radius
    #print "gradient_normal", gradient_normal

    ###Angle of the normal "plane"
    angle_normalplane = numpy.arctan(1/gradient_normal)
    #print "angle_normalplane", angle_normalplane

    ##Angle refraction calculation now
    angle_incidence = incident_rays[:,2] - (numpy.pi / 2) - angle_normalplane
    #print "angle_incidence", angle_incidence
    angle_refracted = numpy.arcsin((spherical_surface[4]*numpy.sin(angle_incidence))/spherical_surface[5])
    #print "angle_refracted", angle_refracted
    #print angle_refracted
    angle_refracted = angle_refracted + (numpy.pi/2.0) + angle_normalplane
    #print "angle_refracted", angle_refracted

    #Outputting final array for refracted rays
    final_ray = numpy.array([intercept_plane_x, intercept_plane_y, angle_refracted]) #the rays are shown vertically
    #Next two lines of code are to make it look the same as the initial incident rays array
    final_ray = numpy.rot90(final_ray, 3)
    refracted_rays = numpy.fliplr(final_ray)

    """
    print "refracted_rays:"
    print refracted_rays
    """

    return refracted_rays

#Task 2
def refraction_2d_det (incident_rays, x_det):
    '''Function to calculate the location of a ray (or multiple rays) as it
    strikes the detector surface.

    Parameters:
        incident_rays: a 2-dimensional numpy array of floats: first index is the
        individual incident light ray, second index for each ray is, the horizontal
        coordinate (metres) of the point of origin of the ray, the corresponding
        vertical coordinate, and the clockwise angle (radians) with respect to the
        positive vertical direction of the direction of the incident ray.

        x_det: a float value specifying the x position of the vertically (y)
        oriented detector surface.

    Output:
        refracted_rays: a 2-dimensional numpy array of floats: individual refracted
        light ray specified by first index, following values specified by second
        index for each ray; the horizontal coordinate (metres) of the point of
        intersection of the corresponding incident ray with the specified detector
        surface, the corresponding vertical coordinate, and the clockwise angle
        (radians).
    '''

    #Creating the final array
    refracted_rays = numpy.zeros(incident_rays.shape, dtype=float)

    #gradients
    #incident ray
    gradient_ray = numpy.tan((numpy.pi/2.0) - incident_rays[:,2])

    #y-intercepts
    #incident ray
    intercept_ray = incident_rays[:,1] - (gradient_ray * incident_rays[:, 0])

    #point of intersection with plane
    #y-coord
    y_intersect_plane = (gradient_ray * x_det) + intercept_ray
    #print y_intersect_plane

    #Outputting final array for refracted rays
    refracted_rays[:,0] = x_det
    refracted_rays[:,1] = y_intersect_plane

    """
    print "refracted_rays:"
    print refracted_rays
    """

    return refracted_rays

#Task 3
def trace_2d (incident_rays, surface_list):
    """Function to calculate the location of a ray (or mutiple rays) as it
    travels through varying surfaces such as planar and spherical, and then
    finally as it hits the detector surface.

    Parameters:
        incident_rays: a 2-dimensional numpy array of floats: first index is the
        individual incident light ray, second index for each ray is, the horizontal
        coordinate (metres) of the point of origin of the ray, the corresponding
        vertical coordinate, and the clockwise angle (radians) with respect to the
        positive vertical direction of the direction of the incident ray.

        surface_list: a list containing the specification of each surface. Each
        element of the list contains nested elements: the first element is a string
        which contains either 'PLA' or 'SPH' or 'DET'. The second element is either
        a float value specifying the x position in the case of DET, or a numpy array
        specifying either a planar or spherical surface in the same formats as used
        for refraction_2d and refration_2d_sph respectively.

    Output:
        refracted_ray_paths: a 3-dimensional array containing the refracted_rays
        output for each surface (float values) in the same order as the surfaces in
        surface_list.
    """

    #Final array storage
    refracted_ray_paths = numpy.zeros((len(surface_list), len(incident_rays), 3), dtype=float)
    #print len(surface_list)
    refracted_ray_paths[0] = incident_rays
    #print refracted_ray_paths

    #number_times was used as a code verifier i.e. is it running?
    number_times = 0

    #For loop which goes through every surface specified by surface_list
    for i in range(len(surface_list)):
        #print i
        #print len(i)
        #surface_type = i[:][0]
        surface_type = surface_list[i][0]
        #print surface_type
        surface_information = surface_list[i][1]
        #print surface_information
        #print i[1]

        #Making sure the code pulls the previous ray
        if i == 0 :
            incident_rays = incident_rays
        else:
            incident_rays = refracted_ray_paths[i-1]

        if surface_type == 'PLA': #Planar surface
            refracted_ray_paths[i] = refraction_2d(incident_rays,surface_list[i][1])
            #print refracted_ray_paths
            #print len(i)
            #print "success PLA"
            number_times = number_times + 1

        if surface_type == 'SPH': #Spherical surface
            refracted_ray_paths[i] = refraction_2d_sph(incident_rays,surface_list[i][1])
            #print refracted_ray_paths
            #print len(i)
            #print "success SPH"
            number_times = number_times + 1

        if surface_type == 'DET': #Detector surface
            refracted_ray_paths[i] = refraction_2d_det(incident_rays,surface_list[i][1])
            #print refracted_ray_paths
            #print len(i)
            #print "success DET"
            number_times = number_times + 1

            break #I'M BREAKING FREE


        #print "number_times", number_times
        #print surface_type
        #print surface_information

    #print refracted_ray_paths

    return refracted_ray_paths

#Task 4
def plot_trace_2d (incident_rays, refracted_ray_paths, surface_list):
    """Function to plot the path of a ray (or multiple rays) as it travels from
    an origin through varying surfaces such as planar and spherical, and then
    finally as it hits the detector surface.

    Parameters:
        incident_rays: a 2-dimensional numpy array of floats: first index is the
        individual incident light ray, second index for each ray is, the horizontal
        coordinate (metres) of the point of origin of the ray, the corresponding
        vertical coordinate, and the clockwise angle (radians) with respect to the
        positive vertical direction of the direction of the incident ray.

        refracted_ray_paths: a 3-dimensional array containing the refracted_rays
        output for each surface (float values) in the same order as the surfaces in
        surface_list.

        surface_list: a list containing the specification of each surface. Each
        element of the list contains nested elements: the first element is a string
        which contains either 'PLA' or 'SPH' or 'DET'. The second element is either
        a float value specifying the x position in the case of DET, or a numpy array
        specifying either a planar or spherical surface in the same formats as used
        for refraction_2d and refration_2d_sph respectively.

    Output:
        a plot containing the surface(s) and ray path(s), a histogram is also
        displayed in the top right, detailing the distribution of the data.

    """

    #Making sure things are defined properly
    refracted_ray_paths = trace_2d(incident_rays,surface_list)
    rrp = refracted_ray_paths
    #print rrp
    surface_list = surface_list

    #Creating the figure to plot the rays and surfaces on
    fig1 = pyplot.figure()

    #Plotting the various surfaces and ray paths
    for i in range(len(surface_list)):
        if surface_list[i][0] == 'SPH': #Plotting spherical surface
            #print surface_list[i][1][6]
            #ax1 = fig1.add_subplot(111, aspect='equal')
            #print surface_list[i][1][0], "test123"

            ss = surface_list[i][1]
            circle_centre_c = ss[2] + (ss[6] / numpy.abs(ss[6])) * numpy.sqrt(((ss[6] ** 2) - (ss[1] ** 2)))
            refraction_2d_sph.ccc = circle_centre_c

            #ax2 = fig1.add_subplot(111, aspect='equal')
            x = numpy.linspace((circle_centre_c - surface_list[i][1][6]),surface_list[i][1][0],100)
            y = numpy.sqrt((surface_list[i][1][6] ** 2) - ((x - circle_centre_c) ** 2))

            pyplot.plot(refraction_2d_sph.ccc, 0)
            #print refraction_2d_sph.ccc
            #print surface_list[i][1][6]
            pyplot.plot(x,(numpy.sqrt((surface_list[i][1][6] ** 2) - ((x - refraction_2d_sph.ccc) ** 2))), 'k-')
            pyplot.plot(x,(-numpy.sqrt((surface_list[i][1][6] ** 2) - ((x - refraction_2d_sph.ccc) ** 2))), 'k-')

        if surface_list[i][0] == 'PLA': #Plotting the planar surfaces
            #fig1 = pyplot.figure()
            #ax3 = fig1.add_subplot(111, aspect='equal')
            pyplot.plot((surface_list[i][1][0], surface_list[i][1][2]), (surface_list[i][1][1], surface_list[i][1][3]), 'k-')

        if surface_list[i][0] == 'DET': #Plotting the detector surfaces
            #print surface_list[i][0]
            #print surface_list[i][1]
            #ax4 = fig1.add_subplot(111, aspect='equal')
            pyplot.axvline(x=surface_list[i][1])

        #Plotting initial segment of the ray path
        if i == 0 and (incident_rays[0][0] < refracted_ray_paths[0][0][0]):
            #print refracted_ray_paths[i][:,0]
            #print refracted_ray_paths[i][:,1]
            #print incident_rays[:,0]
            #print incident_rays[:,1]
            #ax5 = fig1.add_subplot(111, aspect='equal')
            pyplot.plot((incident_rays[:,0], refracted_ray_paths[i][:,0]), (incident_rays[:,1], refracted_ray_paths[i][:,1]), 'k-')

        else: #Plotting the rest of the ray paths
            #ax6 = fig1.add_subplot(111, aspect='equal')
            pyplot.plot((refracted_ray_paths[i-1][:,0], refracted_ray_paths[i][:,0]), (refracted_ray_paths[i-1][:,1], refracted_ray_paths[i][:,1]), 'k-')

    #make code to let x-axis be automatically adjusted according to the detector surface
    #print surface_list[-1]
    det_det = surface_list[-1][1] + 10

    #pyplot.xlim([0,det_det])
    pyplot.xlim(pyplot.xlim()[0], det_det)
    #pyplot.ylim([-10,10]) #defunct as matplotlib auto shapes the figure

    #Histogram plotting
    ax6 = fig1.add_axes([0.72, 0.72, 0.16, 0.16])
    pyplot.hist((refracted_ray_paths[-1][:,1]))
    ax6.axes.get_xaxis().set_visible(True)
    ax6.axes.get_xaxis().set_ticks([])
    ax6.axes.get_xaxis().set_label_text("Distribution")
    ax6.axes.get_yaxis().set_visible(False)

    #ax6.spines['bottom'].set_smart_bounds(False)
    #ax6.xlabel("Histogram")

    #Saving of the diagram
    fig1.savefig('ray_paths.png', dpi=180, bbox_inches='tight')
    #pyplot.autoscale(enable=True, axis='both', tight=None)
    pyplot.show()

#Task 5
def evaluate_trace_2d (refracted_ray_paths, r):
    """Function to evaluate the image quality in the output from trace_2d
    (refracted_ray_paths) on the final detector surface.

    Parameters:
        refracted_ray_paths: a 3-dimensional array containing the refracted_rays
        output for each surface (float values) in the same order as the surfaces in
        surface_list.

        r: a float value specifying the radius of a 1-dimensional 'circle' on the
        detector which is centered around the mean vertical arrival position for
        all the rays.

    Ouput:
        frac: a float value specifying the calculated fraction of rays within the
        specified circle.
    """

    #Finding the average of the y-coords where the rays hit the surface
    ##Refracted ray paths final, the final ray - making sure it is called
    rrpf = refracted_ray_paths
    #Slicing y-coordinates out of final array
    rrpf = refracted_ray_paths[-1][:,1]
    #Calculating the average point
    rrp_mean = numpy.mean(rrpf)

    #Adding r above the average value
    rrp_mean_plus = rrp_mean + r
    #Subtracting r below the average value
    rrp_mean_minus = rrp_mean - r

    #print rrpf, rrp_mean, rrp_mean_plus, rrp_mean_minus

    #Number counts
    ray_circle_yes = 0.0
    ray_circle_no = 0.0

    #Comparing the rays and seeing if they fall within the circle
    for i in rrpf:
        if (i >= rrp_mean_minus) and (i <= rrp_mean_plus):
            ray_circle_yes = ray_circle_yes + 1.0
        else:
            ray_circle_no = ray_circle_no + 1.0
            #the above doesn't really do anything, could be used in the future

    #print "ray_circle_yes", ray_circle_yes

    #Fraction of rays arriving at the final surface
    frac = ray_circle_yes / len(rrpf)

    #print "this is what frac is:",frac

    return frac

#Task 6
def optimize_surf_rad_2d (incident_rays, surface_list, r, n_surf):
    """Function to calculate the optimum radius of curvature values for the
    spherical surfaces which have been specified by the user. The optimum
    radius of curvature values are the ones that give the highest return values
    from evaluate_trace_2d.

    Parameters:
        incident_rays: a 2-dimensional numpy array of floats: first index is the
        individual incident light ray, second index for each ray is, the horizontal
        coordinate (metres) of the point of origin of the ray, the corresponding
        vertical coordinate, and the clockwise angle (radians) with respect to the
        positive vertical direction of the direction of the incident ray.

        surface_list: a list containing the specification of each surface. Each
        element of the list contains nested elements: the first element is a string
        which contains either 'PLA' or 'SPH' or 'DET'. The second element is either
        a float value specifying the x position in the case of DET, or a numpy array
        specifying either a planar or spherical surface in the same formats as used
        for refraction_2d and refration_2d_sph respectively.

        r: a float value specifying the radius of a 1-dimensional 'circle' on the
        detector which is centered around the mean vertical arrival position for
        all the rays.

        n_surf: a 1-dimensional numpy array of integers which specifiy the index
        values of the surfaces in surface_list, the radius of curvature values of
        which is to be optimised.

    Output:
        rad_opt: a 1-dimensional numpy array containing float values of the
        calculated optimal radius values for the specified spherical surfaces
        in n_surf and following the same order.
    """

    rrpr = trace_2d(incident_rays, surface_list)
    pt2d = plot_trace_2d(incident_rays, rrpr, surface_list)
    et2dr = evaluate_trace_2d(rrpr, r)

    #Array to store the optimised radii
    frac_storage = numpy.zeros(1, dtype=float)
    frac_storage[0] = et2dr

    #print rrpr
    #print "evaluated_2d:", et2dr

    #Initial optimisation
    #Creating array to count which radius is which
    rad_opt = numpy.zeros(len(n_surf), dtype=float)
    rad_opt_count = -1.0

    #Array to store the calculated radii
    radius_storage = numpy.zeros(n_surf.shape, dtype=float)

    #Storing the radii that need optimising from the surface list into array
    for k in range(len(n_surf)):
        radius_storage[k] = surface_list[n_surf[k]][1][6]

    #Nested function to help us later in the optimisation
    def function(radius_storage):
        for n in range(len(n_surf)):
            surface_list[n_surf[n]][1][6] = radius_storage[n]
            evaluated_frac = evaluate_trace_2d(trace_2d(incident_rays, surface_list), r)

        return (-1.0 * evaluated_frac)

    frac = evaluate_trace_2d(trace_2d(incident_rays, surface_list), r)

    if frac == 1.0:
        number = 0.0
        for p in n_surf:
            if number >= 0.0 and number <= len(n_surf):
                rad_opt[number] = surface_list[p][1][6]
                number = number + 1.0
            else:
                pass

    else:
        #Nelder-Mead optimisation using scipy
        rad_opt = minimize(function, radius_storage, method="Nelder-Mead", tol=1e-6)
        #print rad_opt

    #return rad_opt.x
    #print rad_opt.x

    #Making sure the radii are returned as wanted
    rad_opt = rad_opt.x

    return rad_opt

if __name__ == '__main__':


    import ray_main
    import numpy

    #Test code the various surfaces, a lot of testing was attempted
    """
    ray_main.refraction_2d(numpy.array([[0.0,5.0, numpy.pi/2.0],[0.0,5.0, 1.1*numpy.pi/2.0]]),
                                        numpy.array([5.0,2.0,6.0,8.0,1.0,1.33]))"""

    """
    ray_main.refraction_2d_sph(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]),
                                            numpy.array([10.0, -10.0,10.0,10.0,1.0,1.33,15.0]))"""

    #ray_main.refraction_2d_det(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]),6.0)

    #ray_main.trace_2d(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]), (["SPH",numpy.array([3.0,4.0,3.0,-4.0,1.4,1.48,5])],["PLA",numpy.array([5.0,2.0,6.0,8.0,1.0,1.33])],["DET",15.0]))

    #a = ray_main.trace_2d(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]), (["SPH",numpy.array([3.0,4.0,3.0,-4.0,1.4,1.48,5])],["PLA",numpy.array([5.0,2.0,6.0,8.0,1.0,1.33])],["DET",15.0]))

    """
    ray_main.trace_2d(numpy.array([[0.0,5.0, numpy.pi/2.0],[0.0,5.0, 1.1*numpy.pi/2.0]]),
                                (["PLA", numpy.array([5.0,2.0,6.0,8.0,1.0,1.33])],
                                ["SPH", numpy.array([6.0,4.0,6.0,-4.0,1.0,1.33,4.0])],
                                ['DET', 15.0]))"""

    """
    ray_main.trace_2d(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]),
                            (["SPH", numpy.array([6.0,4.0,6.0,-4.0,1.0,1.33,4.0])],
                            ["PLA", numpy.array([7.0,2.0,6.0,8.0,1.0,1.33])],
                            ["DET", 6.0]))"""

    """
    print "huh"
    ray_main.plot_trace_2d(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]),
                            trace_2d,
                            (["SPH",numpy.array([4.0,5.0,4.0,-5.0,1.4,1.48,-5])],
                            ["PLA",numpy.array([10.0,3.0,9.0,-7.0,1.0,1.33])],
                            ["DET",16.0]))"""

    #ray_main.plot_trace_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]), trace_2d, ["SPH", numpy.array([0.0, 5.0,0.0,-5.0,1.0,1.33,15.0])])

    """
    a = ray_main.trace_2d(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]),
                        (["SPH",numpy.array([3.0,4.0,3.0,-4.0,1.4,1.48,5])],
                        ["PLA",numpy.array([5.0,2.0,6.0,8.0,1.0,1.33])],
                        ["DET",15.0]))"""

    #ray_main.evaluate_trace_2d(a, 5.0)

    """
    rrp = ray_main.trace_2d(numpy.array([[0.0,1.0,(3.0*numpy.pi)/2.0],[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0],[0.0,4.0,numpy.pi/2]]),[["SPH", numpy.array([20.0, 6.0,20.0,-6.0,1.0,1.33,10.0])], ["SPH", numpy.array([20.0, 6.0,20.0,-6.0,1.33,1.00,-10.0])],
                            ["DET", 30.0]])"""

    """
    #Stupid thin lens
    ray_main.plot_trace_2d(numpy.array([[0.0,1.0,(3.0*numpy.pi)/2.0],[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0],[0.0,4.0,numpy.pi/2]]),rrp,[["SPH", numpy.array([20.0, 6.0,20.0,-6.0,1.0,1.33,10.0])],["SPH", numpy.array([20.0, 6.0,20.0,-6.0,1.33,1.00,-10.0])],
                            ["DET", 30.0]])"""

    """
    #Lucy's test code
    ray_main.plot_trace_2d(numpy.array([[0.0,0.0,numpy.pi/2.0],[0.0,0.0,0.9*numpy.pi/2]]),rrp,[["SPH", numpy.array([10.0, 2.0,10.0,-2.0,1.0,1.52,2.0])],["DET", 16.0]])"""

    """
    ray_main.plot_trace_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0],[0.0,4.0,numpy.pi/2]]),ray_main.trace_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0],[0.0,4.0,numpy.pi/2]]),[["PLA", numpy.array([1.5, -8.0, 2.5, 8.0,1.0,1.33])],["SPH", numpy.array([10.0, 5.0,10.0,-5.0,1.0,1.33,5.0])],["PLA", numpy.array([13.0, -8.0, 15.0, 8.0,1.0,1.33])],["DET", 20.0]]),
                        [["PLA", numpy.array([1.5, -8.0, 2.5, 8.0,1.0,1.33])],["SPH", numpy.array([10.0, 5.0,10.0,-5.0,1.0,1.33,5.0])],["PLA", numpy.array([13.0, -8.0, 15.0, 8.0,1.0,1.33])],["DET", 20.0]])"""

    """
    ray_main.plot_trace_2d(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]),ray_main.trace_2d(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]),(["SPH",numpy.array([10.0,-10.0,10.0,10.0,1.4,1.48, 15.0])],["PLA",numpy.array([5.0,2.0,6.0,8.0,1.0,1.33])],["SPH",numpy.array([12.0,8.0,12.0,-8.0,1.0,1.33,14.0])],["DET",15.0])), (["SPH",numpy.array([10.0,-10.0,10.0,10.0,1.4,1.48, 15.0])],["PLA",numpy.array([5.0,2.0,6.0,8.0,1.0,1.33])],["SPH",numpy.array([12.0,8.0,12.0,-8.0,1.0,1.33,14.0])],["DET",15.0]))"""

    """
    ray_main.plot_trace_2d(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0],[0.0, 2.0, numpy.pi/2.0]]), numpy.array([[[2.0,3.0,1.61],[1.80,2.71,1.75]],[[5.15,2.87,1.64],[5.02,2.13,1.75]],[[6.68,2.77,1.75],[6.31,1.90,1.78]],[[15.0,1.28,0.0],[15.0,0.016,0.0]]]), (["SPH",numpy.array([3.0,4.0,3.0,4.0,1.4,1.48,5])],["PLA",numpy.array([5.0,2.0,6.0,8.0,1.0,1.33])],["SPH",numpy.array([12.0,6.0,12.0,-6.0,1.0,1.33,6.0])],["DET",15.0]))"""

    """
    ray_main.optimize_surf_rad_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0],[0.0,3.0,0.8*numpy.pi/2.0]]),(["SPH",numpy.array([3.0,4.0,3.0,-4.0,1.0,1.33,5])],["SPH",numpy.array([3.0,4.0,3.0,-4.0,1.33,1.0,-5.0])],["DET",15.0]), 1.0, numpy.array([1]))"""

    #ray_main.refraction_2d_sph(numpy.array([[0.0,0.0,numpy.pi/2.0],[0.0,0.0,numpy.pi/4.0]]), numpy.array([10.0, 2.0,10.0,-2.0,1.0,1.52,2.0]))

    """
    ray_main.plot_trace_2d(numpy.array([[0.0,0.0,numpy.pi/2.0],[0.0,0.0,numpy.pi/4.0]]),ray_main.trace_2d(numpy.array([[0.0,0.0,numpy.pi/2.0],[0.0,0.0,numpy.pi/4.0]]), ["SPH", numpy.array([10.0, 2.0,10.0,-2.0,1.0,1.52,2.0])]),[["SPH",numpy.array([10.0, 2.0,10.0,-2.0,1.0,1.52,2.0])], ["DET", 15.0]])"""

    #ray_main.refraction_2d_sph(numpy.array([[0.0,0.0,numpy.pi/2.0],[0.0,0.0,numpy.pi/4.0]]), numpy.array([10.0, 2.0,10.0,-2.0,1.0,1.52,2.0]))

    """
    ray_main.plot_trace_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0],[0.0,4.0,numpy.pi/2]]),ray_main.trace_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0],[0.0,4.0,numpy.pi/2]]),[["PLA", numpy.array([1.5, -8.0, 2.5, 8.0,1.0,1.33])],["SPH", numpy.array([10.0, 5.0,10.0,-5.0,1.0,1.33,5.0])],["PLA", numpy.array([13.0, -8.0, 15.0, 8.0,1.0,1.33])],["DET", 20.0]]),
                        [["PLA", numpy.array([1.5, -8.0, 2.5, 8.0,1.0,1.33])],["SPH", numpy.array([10.0, 5.0,10.0,-5.0,1.0,1.33,5.0])],["PLA", numpy.array([13.0, -8.0, 15.0, 8.0,1.0,1.33])],["DET", 20.0]])"""


    #Rays
    incident_rays = numpy.zeros((1000, 3), dtype=float)
    #print incident_rays
    incident_rays[:,0] = 0.0
    incident_rays[:,1] = numpy.linspace(-3,3,1000)
    incident_rays[:,2] = numpy.pi/2.0

    ray_main.optimize_surf_rad_2d(incident_rays,
                                         [["SPH", numpy.array([30.0, -10.0,30.0,10.0,1.0,1.33,15.0])],["SPH", numpy.array([45.0, 5.0,45.0,-5.0,1.33,1.0,6.0])],["SPH", numpy.array([20.0, 7.0,20.0,-7.0,1.0,2.0,30.0])],
                                          ["DET", 60.0]], 1.0, numpy.array([0,1,2]))
