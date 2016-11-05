"""An optical ray tracing module with a function for dealing with spherical
surfaces as well as planar surfaces, and also with opaque detector surfaces.
Module is also capable of plotting ray paths and surfaces, tgogether with an
illustration of the distribution of rays at the output. Function can be used to
evaluate the image quality of a give optical design, and a function for
optimising the image quality by adjusting radii of survature."""

USER = "Jacky Cao"
USER_ID = "bbvw84"

import numpy
import matplotlib.pyplot as pyplot
import matplotlib.path as pypath
import matplotlib.patches as patches
import matplotlib.hatch as hatch
import sys

#Mini-project original code
def refraction_2d (incident_rays, planar_surface):
    '''A function which can calculate the location of a refracted ray and also
    the angle the refracted ray makes with the vertical. Input requires is an
    array for both incident rays and planar surface.'''

    funcnames = {'ref2d': refraction_2d}

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
    '''[FUCKING EDIT ME] A function which can calculate the location of a refracted ray and also
    the angle the refracted ray makes with the vertical. Input requires is an
    array for both incident rays and planar surface.'''

    #gradients
    #incident ray
    gradient_ray = numpy.tan((numpy.pi / 2.0) - incident_rays[:,2])
    #print "gradient_ray", gradient_ray

    #y-intercepts
    #incident ray
    intercept_ray = (incident_rays[:,1] - (gradient_ray * incident_rays[:, 0]))
    #print "intercept_ray", intercept_ray

    #numpy.seterr(all = 'ignore') #ignores those pesky runtime error
    #critical_angle = numpy.arcsin((spherical_surface[5] / spherical_surface[4]))

    #Critical angle function - corrected
    if spherical_surface[5] >= spherical_surface[4]:
        tir_possible = False
    else:
        tir_possible = True

    #function terminator if the angles given are greater than the critical angle
    if tir_possible and any(angle_incidence) >= critical_angle:
        raise Exception, "at least one incident ray exceeds the critical angle"
        return
    #else:
        #print "Incident angles are valid, function will proceed."

    #Calculating the centre of the circle
    if spherical_surface[6] > 0:
        circle_centre_c = spherical_surface[2] + numpy.sqrt((spherical_surface[6] ** 2) - (spherical_surface[1] ** 2))
        refraction_2d_sph.ccc = circle_centre_c
        ccc = circle_centre_c
        #print "circle_centre_c", circle_centre_c
    else:
        ss = spherical_surface
        circle_centre_c = ss[2] + (ss[6] / numpy.abs(ss[6])) * numpy.sqrt(((ss[6] ** 2) - (ss[1] ** 2)))
        refraction_2d_sph.ccc = circle_centre_c
        ccc = circle_centre_c

    #Working out the point of intersection with surface
    #Quadratic for x coordiante
    ##Determinant calculator
    a = (1 + (gradient_ray ** 2))
    #print "a", a
    b = (2 * intercept_ray * gradient_ray) - (2 * ccc)
    #print "b", b
    c = (ccc ** 2) + (intercept_ray ** 2) - (spherical_surface[6] ** 2)
    #print "c", c
    #Actual Quadratic - using if statement if concace or convex lens
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
    #angle_normalplace = numpy.arctan(- 1 / ((ccc - intercept_plane_x)/(0 - intercept_plane_y)))
    #numpy.arctan((planar_surface[2]-planar_surface[0])/(planar_surface[3]-planar_surface[1]))
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
    '''[FUCKING EDIT ME] A function which can calculate the location of a refracted ray and also
    the angle the refracted ray makes with the vertical. Input requires is an
    array for both incident rays and planar surface.'''

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
    """Work"""

    #Final array storage
    refracted_ray_paths = numpy.zeros((len(surface_list), len(incident_rays), 3), dtype=float)
    #print len(surface_list)
    refracted_ray_paths[0] = incident_rays

    #print refracted_ray_paths

    sl = surface_list
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


        #print "hey look", number_times

        """
        else:
            if surface_type == 'SPH':
                print i[:][0]
                refraction_2d_sph(incident_rays,surface_information)
                print "success2"

                number_times = number_times + 1
        """

        #print surface_type
        #print surface_information
        #print "arse"

        #while surface_type != 'DET':

    #if len(surface_list) > 1:

    #length = len(surface_list)
    #print length

    #print refracted_ray_paths

        #print this_mod(refracted_rays)

    #refracted_ray_paths = numpy.zeros(refracted_rays.shape, dtype=float)

    return refracted_ray_paths

#Task 4
def plot_trace_2d (incident_rays, refracted_ray_paths, surface_list):
    """work"""

    #Making sure things are defined properly
    refracted_ray_paths = trace_2d(incident_rays,surface_list)
    rrp = refracted_ray_paths
    print rrp
    surface_list = surface_list

    fig1 = pyplot.figure()

    #Plotting the various surfaces and ray paths
    for i in range(len(surface_list)):
        if surface_list[i][0] == 'SPH': #Plotting spherical surface
            #print surface_list[i][1][6]
            #ax1 = fig1.add_subplot(111, aspect='equal')
            """
            ax1.add_patch(
                patches.Circle(
                    (refraction_2d_sph.ccc, 0),   # (x,y)
                    surface_list[i][1][6],          # radius
                    fill = False
                )
            )"""

            #print surface_list[i][1][0], "asdjhkasnd"

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

            """
            #Old Cirlce method, secondary method
            ax2.add_patch(
                patches.Arc(
                    xy = [refraction_2d_sph.ccc, 0],
                    width = (2 * surface_list[i][1][6]),
                    height = (2 * surface_list[i][1][6]),
                    angle = 315,
                    theta1 = (numpy.arctan((0 - surface_list[i][1][3])/(refraction_2d_sph.ccc - surface_list[i][1][2])) * (numpy.pi/180)) + 180,
                    theta2 = (numpy.arctan((0 - surface_list[i][1][3])/(refraction_2d_sph.ccc - surface_list[i][1][2])) * (numpy.pi/180)) + 270
                )
            )"""

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
            #print "arse2.0"
            #ax5 = fig1.add_subplot(111, aspect='equal')
            pyplot.plot((incident_rays[:,0], refracted_ray_paths[i][:,0]), (incident_rays[:,1], refracted_ray_paths[i][:,1]), 'k-')

        else: #Plotting the rest of the ray paths
            #ax6 = fig1.add_subplot(111, aspect='equal')
            pyplot.plot((refracted_ray_paths[i-1][:,0], refracted_ray_paths[i][:,0]), (refracted_ray_paths[i-1][:,1], refracted_ray_paths[i][:,1]), 'k-')


    #make code to let x-axis be automatically adjusted according to the detector surface
    #print surface_list[-1]
    det_det = surface_list[-1][1] + 5

    #adjusting y-axis to autmatically size according to the surfaces
    """
    for i in range(len(surface_list)):
        print i
        print surface_list[i][1]

        to_exclude = [-1]
        print surface_list[~numpy.in1d(range(len(surface_list)),to_exclude)]

        max_value = numpy.amax(surface_list[i][1])
        #max_index = surface_list.index(max_value)
        print max_value, "twot" #max_index"""

    #pyplot.xlim([0,det_det])
    pyplot.xlim(pyplot.xlim()[0], det_det)
    #pyplot.ylim([-10,10]) #defunct as matplotlib auto shapes the figure
    fig1.savefig('ray_paths.png', dpi=180, bbox_inches='tight')
    #pyplot.autoscale(enable=True, axis='both', tight=None)
    pyplot.show()

#Task 5
def evaluate_trace_2d (refracted_ray_paths, r):
    """www"""

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

    print rrpf, rrp_mean, rrp_mean_plus, rrp_mean_minus

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

    print "ray_circle_yes", ray_circle_yes

    #Fraction of rays arriving at the final surface
    frac = ray_circle_yes / len(rrpf)

    print frac

    return frac

#Task 6
def optimize_surf_rad_2d (incident_rays, surface_list, r, n_surf):
    """Lol"""

    rrpr = trace_2d(incident_rays, surface_list)
    pt2d = plot_trace_2d(incident_rays, rrpr, surface_list)
    et2dr = evaluate_trace_2d(rrpr, r)



    #print rrpr
    #print et2dr

    for i in n_surf:
        print "WAHAY", i
        print surface_list[i][0]
        print surface_list[i][1][6]
        print ((10.0/100.0) * surface_list[i][1][6])

        #Calculating +ve 10%
        et2dr_plus = surface_list[i][1][6] + ((10.0/100.0) * surface_list[i][1][6])

        #Calculating -ve 10%
        et2dr_minus = surface_list[i][1][6] - ((10.0/100.0) * surface_list[i][1][6])

        print et2dr_plus, et2dr_minus

        if et2dr_plus < et2dr: #or is it a for loop?
            r_new = r + 1.0
            et2dr = evaluate_trace_2d(rrpr, r_new)

            if r_new > 1.0:
                print "no need to further optimise"
                r_accepted = r

            else:
                r_new


    #return rad_opt

if __name__ == '__main__':


    import ray_main
    import numpy

    """
    ray_main.refraction_2d(numpy.array([[0.0,5.0, numpy.pi/2.0],[0.0,5.0, 1.1*numpy.pi/2.0]]),
                                        numpy.array([5.0,2.0,6.0,8.0,1.0,1.33]))"""

    """
    ray_main.refraction_2d_sph(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]),
                                            numpy.array([6.0,4.0,6.0,-4.0,1.0,1.33,4.0]))"""

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
    ray_main.plot_trace_2d(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]),
                            trace_2d,
                            (["SPH",numpy.array([3.0,4.0,3.0,-4.0,1.4,1.48,5])],
                            ["PLA",numpy.array([5.0,2.0,6.0,8.0,1.0,1.33])],
                            ["DET",15.0]))"""

    #ray_main.plot_trace_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0]]), trace_2d, ["SPH", numpy.array([0.0, 5.0,0.0,-5.0,1.0,1.33,15.0])])

    """
    a = ray_main.trace_2d(numpy.array([[0.0,3.0, numpy.pi/2.0],[0.0,3.0, 1.1*numpy.pi/2.0]]),
                        (["SPH",numpy.array([3.0,4.0,3.0,-4.0,1.4,1.48,5])],
                        ["PLA",numpy.array([5.0,2.0,6.0,8.0,1.0,1.33])],
                        ["DET",15.0]))"""

    #ray_main.evaluate_trace_2d(a, 5.0)

    ray_main.optimize_surf_rad_2d(numpy.array([[0.0,5.0,numpy.pi/2.0],[0.0,5.0,1.1*numpy.pi/2.0],[0.0,3.0,0.8*numpy.pi/2.0]]),
                                     [["PLA", numpy.array([4.5, -8.0, 5.5, 8.0,1.0,1.33])],
                                      ["SPH", numpy.array([10.0, -10.0,10.0,10.0,1.0,1.33,15.0])],
                                      ["DET", 15.0]], 1.0, numpy.array([1]))
