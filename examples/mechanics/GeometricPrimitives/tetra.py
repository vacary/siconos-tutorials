#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import siconos.numerics as sn
import siconos.kernel as sk

from siconos.mechanics.collision.convexhull import ConvexHull


# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a tetrahedron as a convex shape
    import numpy
    pts = numpy.array([(-1.0, 1.0, 1.0),
                       (1.0, -1.0, 1.0),
                       (-1.0, -1.0, 1.0),
                       (0.0, 0.0, -1.0)])
    io.add_convex_shape('Tetra', pts - pts.mean(0),
                        insideMargin=0.01, outsideMargin=0.0)

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (10, 10, 1),
                           insideMargin=0.01, outsideMargin=0.0)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.01, e=0.7,
                                      collision_group1=1,
                                      collision_group2=2)

    # computation of inertia and volume
    ch = ConvexHull(pts)
    inertia, volume = ch.inertia(ch.centroid())

    # The tetra object made with an unique Contactor : the tetrahedron
    # shape.  As a mass is given, it is a dynamic system involved in
    # contact detection and in the simulation.  With no group id
    # specified the Contactor belongs to group 0
    io.add_object('tetra', [Contactor('Tetra', collision_group=1)],
                  translation=[0, 0, 4],
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground', collision_group=2)],
                  translation=[0, 0, -0.1])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

from siconos.mechanics.collision.bullet import SiconosBulletOptions

# Control the number of perturbations applied to generate multipoint
# surface-surface contact manifolds.  Default is 5 and 5, this is
# just demonstrating how to change them.
bullet_options = SiconosBulletOptions()
bullet_options.perturbationIterations = 4
bullet_options.minimumPointsPerturbationThreshold = 4

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8


with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.

    # print(pydoc.render_doc(io.run, "Help on %s"))

    io.run(with_timer=False,
           bullet_options=bullet_options,
           t0=0,
           T=20,
           h=0.005,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=None)
