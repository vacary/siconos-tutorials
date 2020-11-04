#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import siconos.numerics as sn
import siconos.kernel as sk
# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    R = 0.2
    # Definition of a sphere
    io.add_primitive_shape('Sphere', 'Sphere', (R,),
                           insideMargin=0.2, outsideMargin=0.0)

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (10, 10, 0.1),
                           insideMargin=0.05, outsideMargin=0.0)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_binary_cohesive_nsl('contact', mu=1.0, e=0.0, sigma_c=1e+01, delta_c=1e-01)

    n_spheres=2
    
    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    for i in range(n_spheres):
        for j in range(n_spheres -i):
            io.add_object('sphere' + str(i)+str(j), [Contactor('Sphere')],
                          translation=[0, j*2*R, 2*R+i*2*R],
                          velocity=[0, 0, 0, 0, 0, 0],
                          mass=1)
        for j in range(1,n_spheres -i):
            io.add_object('sphere' + str(i)+str(-j), [Contactor('Sphere')],
                          translation=[0, -j*2*R, 2*R+i*2*R],
                          velocity=[0, 0, 0, 0, 0, 0],
                          mass=1)
    for i in range(1,n_spheres):
        for j in range(n_spheres -i):
            io.add_object('sphere' + str(-i)+str(j), [Contactor('Sphere')],
                          translation=[0, j*2*R, 2*R-i*2*R],
                          velocity=[0, 0, 0, 0, 0, 0],
                          mass=1)
            for j in range(1,n_spheres -i):
                io.add_object('sphere' + str(-i)+str(-j), [Contactor('Sphere')],
                              translation=[0, -j*2*R, 2*R-i*2*R],
                              velocity=[0, 0, 0, 0, 0, 0],
                              mass=1)
    # for i in range(n_spheres):
    #     io.add_object('sphere' + str(-i), [Contactor('Sphere')],
    #                   translation=[0, 0 , 2*R+i*2*R],
    #                   velocity=[0, 0, 0, 0, 0, 0],
    #                   mass=1)
 
    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 0, -2*R*n_spheres-0.05])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

from siconos.mechanics.collision.bullet import SiconosBulletOptions
bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.04

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-4

h=1e-3
T=1000*h
with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(with_timer=True,
           bullet_options=bullet_options,
           face_class=None,
           edge_class=None,
           t0=0,
           T=T,
           h=h,
           multipoints_iterations=False,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=1)
