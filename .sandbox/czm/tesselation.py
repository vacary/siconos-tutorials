#!/usr/bin/env python

from siconos.io.mechanics_run import MechanicsHdf5Runner, MechanicsHdf5Runner_run_options
import siconos.numerics as sn
import siconos.kernel as sk

import read_tess

import sys
import numpy

if (len(sys.argv) < 2):
    dist = 'uniform'
    mu = 0.1
else:
    dist = sys.argv[1]
    mu = sys.argv[2]

if not dist in ['uniform', 'double', 'exp']:
    print("dist = [uniform | double | exp]")
    sys.exit(1)
if float(mu) < 0.1 or float(mu) > 2.0:
    print("mu = [0.1 .. 2.0]")
    sys.exit(1)


fn = 'tesselation-{0}-mu-{1}.hdf5'.format(dist,mu)


test=False

if test:
    hstep=1e-3
    T = hstep*1000
    filename='n5-id1.tess'
    fn = 'tesselation.hdf5'
else:
    hstep=1e-3
    T = hstep*10000
    filename='n100-id1.tess'
    fn = 'tesselation.hdf5'

tesselation = read_tess.read_tesselation(filename)

#print(tesselation)
#input()
from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull

with MechanicsHdf5Runner(mode='w', io_filename=fn) as io:

    all_faces = tesselation['face']
    all_edges = tesselation['edge']
    all_vertices = tesselation['vertex']
    #print(all_faces[:])
    #input()
    for p in tesselation['polyhedron']:
        #print('p=', p)
        #create set of vertices
        faces= p[1]
        vertices=[]
        for f in faces:
            #print('f', f)
            [face]  = [x for x in all_faces if x[0] == abs(f)]
            #print(face)
            vertice = face[1]
            vertices.extend(vertice)
        vertices=list(set(vertices))
        #print('ver', vertices)
        if (len(vertices) <=  1):
            print('the number of vertices must be mroe than 0')
        vertices_coordinates=[]
        for v in vertices:
            #print(v)
            [vertex]  = [x for x in all_vertices if x[0] == abs(v)]
            #print(vertex)
            vertices_coordinates.append([vertex[1],vertex[2],vertex[3]])
        #print(vertices_coordinates)

        # Definition of a polyhedron as a convex shape
        cname = 'polyhedron' + str(p)
        ch = ConvexHull(vertices_coordinates)
        cm = ch.centroid()
        vertices_coordinates = (numpy.array(vertices_coordinates)[:] - cm[:])
        io.add_convex_shape(cname, vertices_coordinates, insideMargin=0.0)

        ch = ConvexHull(vertices_coordinates)
        inertia,volume=ch.inertia(ch.centroid())

        density=2300
        #print('geometric inertia:', inertia)
        #print('volume:', volume)
        # print('mass:', volume*density)
        # print('inertia:', inertia*density)


        name = 'polyhedron_bdy' + str(p)

        io.add_object(name,
                      [Contactor(cname, collision_group=1)],
                      translation=cm*1.0,
                      orientation=[1.0, 0.0, 0., 0.],
                      velocity=[0.,0.,0.,0.,0.,0.],
                      mass=volume*density,
                      inertia=inertia*density)
        #input()

    import math
    angle = math.pi/4.0

    
    
        
    io.add_primitive_shape('Ground', 'Box', (4.0, 8.0, 0.2),
                           insideMargin=0.0, outsideMargin=0.0)
    io.add_object('ground', [Contactor('Ground', collision_group=0)],
                  translation=[.0, .0, -0.5-0.5],
                  orientation = [math.cos(angle/2), math.sin(angle/2), math.sin(angle/2), math.sin(angle/2) ])


    #io.add_Newton_impact_friction_nsl('contact_brick_brick', mu=1.0, e=0.0, collision_group1=1, collision_group2=1)
    #io.add_Newton_impact_friction_nsl('contact_ground_brick', mu=1.0, e=0.0, collision_group1=0, collision_group2=1)
    
    io.add_binary_cohesive_nsl('contact_brick_brick', mu=1.0, e=0.0, sigma_c=1.5e+03, delta_c=1e-02, collision_group1=1, collision_group2=1 )
    io.add_Newton_impact_friction_nsl('contact_ground_brick', e=0.0, mu=1.0, collision_group1=0, collision_group2=1)



# Create solver options
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3


run_options=MechanicsHdf5Runner_run_options()
run_options['t0']=0
run_options['T']=T
run_options['h']=hstep
run_options['theta']= 1.0

#run_options['bullet_options']=bullet_options
run_options['solver_options']=options
run_options['constraint_activation_threshold']=1e-04

run_options['Newton_options']=sk.SICONOS_TS_LINEAR
#run_options['skip_last_update_output']=True
#run_options['skip_reset_lambdas']=True

run_options['osns_assembly_type']= sk.REDUCED_DIRECT

run_options['verbose']=True
run_options['with_timer']=True
run_options['explode_Newton_solve']=True
run_options['explode_computeOneStep']=True

#run_options['numerics_verbose']=True
#run_options['numerics_verbose_level']=0

run_options['output_frequency']=1

#run_options['time_stepping']=None
#run_options['Newton_max_iter']=1
#run_options['output_contact_index_set']=1



with MechanicsHdf5Runner(mode='r+', io_filename=fn) as io:
    io.run(run_options)
