# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
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
import numpy as np
import os
import sys

# siconos
from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull

# from convexhull_modif import ConvexHull
from siconos.mechanics.collision.bullet import SiconosBulletOptions
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
)
import siconos.numerics as sn
import siconos.integrators

# stl mesh
# from stl import mesh
import trimesh

# function generate a shape
import generate_shape as gen_shape

import random
import argparse

random.seed(42)
###############################################################################
# Definition
###############################################################################

# Script name
sName = os.path.basename(sys.argv[0])[:-3]

# Working DIrertory
workDir = os.path.dirname(os.path.realpath(sys.argv[0]))

# Simu parameters---------------------------------------------------------------

# time
time = 20.0  # sec  ok avec 1 dt


# Soil param -------------------------------------------------------------------
# concat

zone_id = [1, 2]

nb_zones = len(zone_id)

# zone 1:
e1 = 0.0
mu1 = 0.8
mur1 = 0.4

# zone 2:
e2 = 0.0
mu2 = 0.7
mur2 = 0.1

e_c = [e1, e2]
mu_c = [mu1, mu2]
mu_r_c = [mu1, mur2]

# blocks
density = 2600
# volumes
Vblock_min = 1.0
Vblock_max = 2.0
# shapes
elBlocx = 1.1
elBlocy = 1.2
nbPts = 40

# height_fall
height_fall_min = 1.0
height_fall_max = 2.0
mnt_file_prefix = "./data/mnt_nord_sud_id"
raster_dep_file = "./data/zone_dep_mnt_nord_sud.stl"
mnt_raster_file = "./data/mnt_nord_sud.asc"

# passage argument --output

parser = argparse.ArgumentParser(description="Simulation Rockfall")

nb_rocks_default = 10
parser.add_argument(
    "--nblocks",
    type=int,
    default=nb_rocks_default,
    help="number of blocks for each run",
)

parser.add_argument(
    "--output",
    type=str,
    default="simu_rockfall_biplan.hdf5",
    help="Full path to the result file (hdf5)",
)
args = parser.parse_args()
# nom de fichier hdf5
fn = args.output
nb_rocks = args.nblocks
time += nb_rocks_default

# time = 14.2
# output = f"time is {time}..."
# input(f"time is {time}...")


class deleterock:
    def __init__(self):
        pass

    def initialize(self, io):
        self._io = io
        pass

    def call(self, step):
        print("call death hook at step", step)
        positions = self._io._io.positions(self._io._nsds)
        velocities = self._io._io.velocities(self._io._nsds)
        # list contact points
        contact_points = self._io._io.contactPoints(self._io._nsds, 1)
        if contact_points is not None:

            for cp in contact_points:
                # print("cp", cp)
                inter_id = int(cp[22])
                print("inter_id", inter_id)
                ds1_id = int(cp[23])
                ds2_id = int(cp[24])
                print("ds id :", ds1_id, ds2_id)
                contact_normal = np.array([cp[7], cp[8], cp[9]])
                contact_force = np.array([cp[10], cp[11], cp[12]])
                fn = np.dot(contact_normal, contact_force)
                ft = np.linalg.norm(
                    contact_force
                    - np.dot(contact_normal, contact_force) * contact_normal
                )
                print("positions", positions)
                print("velocities", velocities)

                print("contact_normal", contact_normal)
                print("contact_force", contact_force)
                print("fn", fn)
                print("ft", ft)

                if np.linalg.norm(fn) > 0:
                    print("ft/fn", ft / fn)
                if np.isnan(velocities).any():
                    print("velocities", np.array(velocities)[:, ds1_id])
                    input()

                """
                if ds1_id == ds2_id:
                    print("interaction with a static object")
                    #inter = self._io._nsds.interaction(inter_id)
                    inter = siconos.modeling.interactions
                    contact_r = inter.relation()
                    print("contact_r.bodyShapeRecordA", contact_r.distance())
                    print("contact_r.bodyShapeRecordA", contact_r.bodyShapeRecordA)
                    print("contact_r.bodyShapeRecordB", contact_r.bodyShapeRecordB)

                    print("contact_r.bodyShapeRecordB.staticBody")
                    print(
                        "contact_r.bodyShapeRecordB.staticBody",
                        contact_r.bodyShapeRecordB.staticBody,
                    )
                    print(
                        "contact_r.bodyShapeRecordB.staticBody.number",
                        contact_r.bodyShapeRecordB.staticBody.number,
                    )
                    print(
                        "contact_r.bodyShapeRecordB.ds", contact_r.bodyShapeRecordB.ds
                    )

                    print("contact_r.bodyShapeRecordA.staticBody")
                    print(
                        "contact_r.bodyShapeRecordA.staticBody",
                        contact_r.bodyShapeRecordA.staticBody,
                    )
                    print(
                        "contact_r.bodyShapeRecordA.staticBody.ds.number()",
                        contact_r.bodyShapeRecordA.ds.number(),
                    )

                    if contact_r.bodyShapeRecordB.staticBody.number == static_body_id:
                        self._io._interman.removeStaticBody(
                            contact_r.bodyShapeRecordB.staticBody
                        )
                        # remove the body from that list of static object completely
                        for s in self._io._static:
                            # print(self._io._static[s])
                            # print(self._io._static[s]['number'])
                            if self._io._static[s]["number"] == static_body_id:
                                s_remove = s
                        self._io._static.pop(s_remove)
        else:
            # print('no contact points')
            pass

            # second way (faster) :  direct access to nsds positions
        # positions  = self._io._io.positions(self._io._nsds)
        # if positions is not None:
        #     z = positions[:,3]
        #     # We search for the ds index that are below a given criteria
        #     ds_idx = numpy.nonzero(z < -2)[0]
        #     for i in ds_idx :
        #         n_ds = int(positions[i,0])
        #         ds = self._io._nsds.dynamicalSystem(n_ds)
        #         self._io._interman.removeBody(ds)
        #         self._io._nsds.removeDynamicalSystem(ds)
        #         print('remove ds number ', ds.number(), ' with height = ', z[i])
        """


###############################################################################
# Simulation
###############################################################################

with MechanicsHdf5Runner(io_filename=fn) as io:

    # soil ---------------------------------------------------------------------
    print("Resample Soil mask... \n")
    for i in range(0, nb_zones):
        dem_file = mnt_file_prefix + str(zone_id[i]) + ".stl"
        # if i==1:
        # dem_file='mnt_2m_PR_simple_'+str(i)+'_mod.stl'
        # else:
        # dem_file='mnt_2m_PR_simple_'+str(i)+'.stl'
        mesh1 = io.add_mesh_from_file(
            "dem_TIN_shape%d" % i,
            dem_file,
            scale=1,
            insideMargin=0.0,
            outsideMargin=0.0,
        )
        io.add_object(
            "dem_TIN%d" % i,
            [Contactor("dem_TIN_shape%d" % i, collision_group=i)],
            translation=[0.0, 0.0, 0.0],
        )

    # blocks

    # chargement zones depart
    mesh = trimesh.load(raster_dep_file)
    areas = mesh.area_faces
    prob = areas / areas.sum()

    # creation blocks
    for i in range(0, nb_rocks):
        print("n_rock", i)
        nameDs = "block" + str(i)
        nameShape = "block" + str(i) + "-shape"
        initial_angles = np.random.uniform(0, 2 * np.pi, 3)
        (c1, c2, c3), (s1, s2, s3) = np.cos(initial_angles), np.sin(initial_angles)
        initial_orientation = [
            c1 * c2 * c3 - s1 * s2 * s3,
            s1 * s2 * c3 + c1 * c2 * s3,
            s1 * c2 * c3 + c1 * s2 * s3,
            c1 * s2 * c3 - s1 * c2 * s3,
        ]

        # Position initiale
        height_fall = height_fall_min + (
            height_fall_max - height_fall_min
        ) ** np.random.uniform(0, 1, 1)
        triangle_index = np.random.choice(len(mesh.faces), p=prob)
        triangle = mesh.triangles[triangle_index]
        A, B, C = triangle
        r1 = np.sqrt(np.random.rand())
        r2 = np.random.rand()
        point = (1 - r1) * A + r1 * (1 - r2) * B + r1 * r2 * C
        Xpos = point[0]
        Ypos = point[1]
        Zpos = point[2] + Vblock_max**0.33 + height_fall[0]

        # posInit=[-0.5*heightmap.shape[0]+indx_dep[0][0],-0.5*heightmap.shape[1]+indx_dep[0][1],315]
        posInit = [Xpos, Ypos, Zpos]
        dest_vol = Vblock_min + np.random.rand() * (Vblock_max - Vblock_min)
        vertices = gen_shape.generate_shape(nbPts, elBlocx, elBlocy, dest_vol)
        """
        vertices = np.array([
            [-1, -1, -1],
            [-1, -1,  1],
            [-1,  1, -1],
            [-1,  1,  1],
            [ 1, -1, -1],
            [ 1, -1,  1],
            [ 1,  1, -1],
            [ 1,  1,  1]
        ], dtype=float)
        """
        # create the convexhull

        ch = ConvexHull(vertices)
        cm = ch.centroid()
        # move the vertices to center the center of mass at 0.0
        vertices = np.array(vertices)[:] - cm[:]
        # ch = ConvexHull(vertices)
        # cm = ch.centroid()
        inertia, area = ch.inertia(cm)
        mass_block = 2650.0 * 2.0

        inertia = inertia * mass_block
        # random block shape
        io.add_convex_shape(nameShape, vertices, outsideMargin=0.0)
        io.add_object(
            nameDs,
            [Contactor(nameShape, collision_group=100)],
            translation=posInit,
            velocity=[0, 0, 0, 0, 0, 0],
            orientation=initial_orientation,
            mass=mass_block,
            inertia=inertia,
            time_of_birth=float(i),
            time_of_death=time + float(i),
        )

    # contact laws
    for i in range(0, nb_zones):
        io.add_Newton_impact_rolling_friction_nsl(
            "contact_soil_%d" % int(i),
            e=e_c[i - 1],
            mu=mu_c[i - 1],
            mu_r=mu_r_c[i - 1],
            collision_group1=100,
            collision_group2=i,
        )

# TEST OPTIONS DE BASE EXAMPLE TUTORIALS

np_timestep = 0.001


bullet_options = SiconosBulletOptions()
# bullet_options.worldScale = .1
# bullet_options.contactBreakingThreshold = 0.4


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sn.solver_options_create(sn.solver_ids.SICONOS_ROLLING_FRICTION_3D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-5


run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = time
run_options["h"] = np_timestep

run_options["bullet_options"] = bullet_options
run_options["solver_options"] = options


# run_options['skip_last_update_output']=True
# run_options['skip_reset_lambdas']=True
# run_options['osns_assembly_type']= siconos.nonsmooth_formulations.REDUCED_DIRECT

run_options["Newton_options"] = siconos.simulation.NONLINEAR

run_options["Newton_max_iter"] = 50
run_options["Newton_tolerance"] = 1e-8

run_options["verbose"] = True
run_options["with_timer"] = False
# run_options["explode_computeOneStep_in_python"] = True
run_options["explode_computeOneStepNSProblem_in_python"] = True

# run_options['violation_verbose'] = True
run_options["output_frequency"] = 1

run_options["time_stepping"] = None

dr = deleterock()

run_options["end_run_iteration_hook"] = dr

with MechanicsHdf5Runner(mode="r+", io_filename=fn) as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.

    # We run the simulation with relatively low precision because we
    # are interested mostly just in the location of contacts with the
    # height field, but we are not evaluating the performance of
    # individual contacts here.
    io.run(run_options)
