#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import siconos.numerics as sn
import siconos.kernel as sk

import math

def diamond_sphere_pack(io,R, n_spheres):
    initial_gap=2*R*n_spheres
    for i in range(n_spheres):
        for j in range(n_spheres -i):
            io.add_object('sphere' + str(i)+str(j), [Contactor('Sphere')],
                          translation=[0, j*2*R, initial_gap+i*2*R],
                          velocity=[0, 0, 0, 0, 0, 0],
                          mass=1)
        for j in range(1,n_spheres -i):
            io.add_object('sphere' + str(i)+str(-j), [Contactor('Sphere')],
                          translation=[0, -j*2*R, initial_gap+i*2*R],
                          velocity=[0, 0, 0, 0, 0, 0],
                          mass=1)
    for i in range(1,n_spheres):
        for j in range(n_spheres -i):
            io.add_object('sphere' + str(-i)+str(j), [Contactor('Sphere')],
                          translation=[0, j*2*R, initial_gap-i*2*R],
                          velocity=[0, 0, 0, 0, 0, 0],
                          mass=1)
            for j in range(1,n_spheres -i):
                io.add_object('sphere' + str(-i)+str(-j), [Contactor('Sphere')],
                              translation=[0, -j*2*r, initial_gap-i*2*R],
                              velocity=[0, 0, 0, 0, 0, 0],
                              mass=1)
                
def square_sphere_pack(io,R, n_spheres, rotation_angle=0):
    initial_gap=R/5.
    initial_gap=0.0
    for i in range(n_spheres):
        for j in range(n_spheres):
            center_y = j*2*R
            center_z = initial_gap+i*2*R+R

            y = math.cos(rotation_angle) * center_y  -math.sin(rotation_angle)* center_z
            z = math.sin(rotation_angle) * center_y  +math.cos(rotation_angle)* center_z
            

            
            io.add_object('sphere' + str(i)+str(j), [Contactor('Sphere',collision_group=1)],
                          translation=[0, y, z],
                          velocity=[0, 0, 0, 0, 0, 0],
                          mass=1)
            
def cube_sphere_pack(io,R, n_spheres, rotation_angle=0):
    initial_gap=R/5.
    initial_gap=0.0
    for i in range(n_spheres):
        for j in range(n_spheres):
            for k in range(n_spheres):
                center_x = i*2*R
                center_y = j*2*R
                center_z = initial_gap+k*2*R+R
                x = center_x 
                y = math.cos(rotation_angle) * center_y  -math.sin(rotation_angle)* center_z
                z = math.sin(rotation_angle) * center_y  +math.cos(rotation_angle)* center_z
            

            
                io.add_object('sphere' + str(i)+str(j)+str(k), [Contactor('Sphere',collision_group=1)],
                              translation=[x, y, z],
                              velocity=[0, 0, 0, 0, 0, 0],
                              mass=1)
   



# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    R = 0.2
    # Definition of a sphere
    io.add_primitive_shape('Sphere', 'Sphere', (R,),
                           insideMargin=0.2, outsideMargin=0.0)
    R_impact = 1.0
    io.add_primitive_shape('Sphere_impact', 'Sphere', (R_impact,),
                           insideMargin=0.2, outsideMargin=0.0)


    ground_thickness=0.2
    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (10, 10, ground_thickness),
                           insideMargin=0.05, outsideMargin=0.0)


    n_spheres=6
    angle=math.pi/6
    angle=0.0
    cube_sphere_pack(io, R, n_spheres,rotation_angle=angle)

    # for i in range(n_spheres):
    #     io.add_object('sphere' + str(-i), [Contactor('Sphere')],
    #                   translation=[0, 0 , 2*R+i*2*R],
    #                   velocity=[0, 0, 0, 0, 0, 0],
    #                   mass=1)

    # io.add_object('sphere_impact' , [Contactor('Sphere_impact',collision_group=0)],
    #               translation=[n_spheres*R-R, n_spheres*R-R, n_spheres*2*R+2*R+R_impact],
    #               velocity=[0, 0, -1, 0, 0, 0],
    #               mass=1)
 

    
    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground', collision_group=0)],
                  translation=[0, 0, -ground_thickness/2.0])
    
    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_binary_cohesive_nsl('contact_sphere_sphere', mu=0.0, e=0.0, sigma_c=1.5e-01, delta_c=1e-02, collision_group1=1, collision_group2=1 )
    
    #io.add_Newton_impact_friction_nsl('contact_sphere_sphere', e=0.0, mu=1.0, collision_group1=1, collision_group2=1)
    
    io.add_Newton_impact_friction_nsl('contact_ground_sphere', e=0.0, mu=1.0, collision_group1=0, collision_group2=1)

    

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

from siconos.mechanics.collision.bullet import SiconosBulletOptions
bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 0.04

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3
options.iparam[sn.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 50
internal_options = sk.solver_options_get_internal_solver(options, 0)
internal_options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 10


h=1e-3

T=347*h
T=4000*h
#T=10*h
#T=1000*h
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
           output_frequency=10)
