# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2024 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
# spheres in a box of size:
# mkspheres.lx*mkspheres.ly*mspheres.lz
# for n spheres:
# ./mkspheres.py <n>
# ./spheres_in_a_box.py coors-<n>.txt radii-<n>.txt
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
from siconos.mechanics.collision.tools import Contactor

from siconos.io.FrictionContactTrace import FrictionContactTraceParams

import siconos.numerics as sn

from siconos.mechanics.collision.bullet import SiconosBulletOptions
from math import pi
import numpy
import sys
import mkspheres

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1000
bullet_options.contactBreakingThreshold = 0.0002
bullet_options.perturbationIterations = 0
bullet_options.minimumPointsPerturbationThreshold = 0

hstep = 0.001
theta = 0.50001
itermax = 1000
tolerance = 1e-7
lx = mkspheres.lx
ly = mkspheres.ly
lz = mkspheres.lz
margin_ratio = 1.0e-5
wthick = mkspheres.lz / 10
zoffset = -wthick / 2

margin_max = margin_ratio * mkspheres.radius_max
if len(sys.argv) > 1:
    coors_filename = sys.argv[1]
    radii_filename = sys.argv[2]
else:
    coors_filename = "coors-18.txt"
    radii_filename = "radii-18.txt"

coors = numpy.loadtxt(coors_filename)
radii = numpy.loadtxt(radii_filename)

if len(radii.shape) == 0:
    radii = [radii]

nb_laid_particles = len(radii)

print(nb_laid_particles)

solver = sn.solver_ids.SICONOS_FRICTION_3D_NSGS
fileName = "spheres-in-a-box-{0}".format(nb_laid_particles)
title = "SpheresBox"
description = """
Spheres in a box, generation with lmgc90 granulo_Random
number of spheres: {0}
radius min       : {1}
radius max       : {2}
box size x       : {3}
box size y       : {4}
box size z       : {5}
Moreau TimeStepping: h={6}, theta={7}
One Step non smooth problem: {8}, maxiter={9}, tol={10}
""".format(
    nb_laid_particles,
    mkspheres.radius_min,
    mkspheres.radius_max,
    lx,
    ly,
    lz,
    hstep,
    theta,
    solver,
    itermax,
    tolerance,
)

mathInfo = ""

friction_contact_trace_params = FrictionContactTraceParams(
    dump_itermax=10000,
    dump_probability=None,
    fileName=fileName,
    title=title,
    description=description,
    mathInfo=mathInfo,
)

# Create solver options
options = sn.solver_options_create(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = itermax
options.dparam[sn.params.SICONOS_DPARAM_TOL] = tolerance

with MechanicsHdf5Runner(io_filename="siab-{0}.hdf5".format(nb_laid_particles)) as io:
    # Definition of the ground shape
    io.add_primitive_shape(
        "Ground",
        "Box",
        (lx, ly, wthick),
        insideMargin=margin_max,
        outsideMargin=margin_max,
    )
    io.add_primitive_shape(
        "Wall",
        "Box",
        (lx, ly, wthick),
        insideMargin=margin_max,
        outsideMargin=margin_max,
    )

    io.add_object("ground", [Contactor("Ground")], translation=[0, 0, zoffset])
    io.add_object(
        "wall1",
        [Contactor("Wall")],
        translation=[-lx / 2 - wthick / 2, 0, lz / 2 + zoffset],
        orientation=([0, 1, 0], pi / 2.0),
    )
    io.add_object(
        "wall2",
        [Contactor("Wall")],
        translation=[lx / 2 + wthick / 2, 0, lz / 2 + zoffset],
        orientation=([0, 1, 0], pi / 2.0),
    )
    io.add_object(
        "wall3",
        [Contactor("Wall")],
        translation=[0, -ly / 2 - wthick / 2, lz / 2 + zoffset],
        orientation=([1, 0, 0], pi / 2.0),
    )
    io.add_object(
        "wall4",
        [Contactor("Wall")],
        translation=[0, ly / 2 + wthick / 2, lz / 2 + zoffset],
        orientation=([1, 0, 0], pi / 2.0),
    )

    io.add_Newton_impact_friction_nsl("contact", mu=0.1, e=0.0)
    for i in range(nb_laid_particles):
        rad = radii[i]
        margin = rad * margin_ratio
        io.add_primitive_shape(
            "Sphere-{0}".format(i),
            "Sphere",
            (rad,),
            insideMargin=margin,
            outsideMargin=margin,
        )
        mass = (4.0 / 3) * pi * (radii[i] ** 3) * 2320
        I_ = (2.0 / 5) * mass * (radii[i] * radii[i])
        io.add_object(
            "sphere-{0}".format(i),
            [Contactor("Sphere-{0}".format(i))],
            translation=[coors[3 * i], coors[3 * i + 1], coors[3 * i + 2]],
            velocity=[0, 0, 0, 0, 0, 0],
            inertia=[I_, I_, I_],
            mass=mass,
        )

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = 8.
run_options["h"] = hstep


run_options["solver_options"] = options
run_options["bullet_options"] = bullet_options
# run_options['constraint_activation_threshold']=1e-05


run_options["Newton_max_iter"] = 1
run_options["output_frequency"] = None

# run_options["verbose"] = False
run_options["with_timer"] = False
# run_options["violation_verbose"] = True

run_options['numerics_verbose'] = False
run_options['numerics_verbose_level'] = 0

run_options["with_timer"]=True
run_options["explode_computeOneStep_in_python"]=True

run_options["friction_contact_trace_params"]=friction_contact_trace_params

with MechanicsHdf5Runner(
    io_filename="siab-{0}.hdf5".format(nb_laid_particles), mode="r+"
) as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(run_options)
