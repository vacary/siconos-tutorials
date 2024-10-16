#!/usr/bin/env python

#
# Example of a bouncing ball with OpenCascade contactors and occ distance
#

from siconos.mechanics.collision.tools import Volume, Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
from siconos import numerics
import siconos.numerics as sn
import siconos.kernel as sk
import siconos.io.mechanics_run
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCC.gp import gp_Pnt



siconos.io.mechanics_run.set_backend('occ')

sphere = BRepPrimAPI_MakeSphere(1.).Shape()
ground = BRepPrimAPI_MakeBox(gp_Pnt(-50, -50, 0), 100., 100., .5).Shape()

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    io.add_occ_shape('Sphere', sphere)
    io.add_occ_shape('Ground', ground)

    io.add_object('sphere',
                  [Volume('Sphere'),
                   Contactor('Sphere', contact_type='Face', contact_index=0)],
                  mass=1, translation=[0, 0, 10], velocity=[0, 0, 0, 0, 0, 0])

    io.add_object('ground',
                  [Contactor('Ground', contact_type='Face', contact_index=5)],
                  translation=[0, 0, 0])

    io.add_interaction('sphere-ground',
                       'sphere', 'Sphere-1',
                       'ground', 'Ground-0',
                       distance_calculator='occ',
                       offset1=0.01)

    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0.9)
    
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.

    io.run(with_timer=False,
           gravity_scale=1,
           t0=0,
           T=10,
           h=0.0005,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=None)
