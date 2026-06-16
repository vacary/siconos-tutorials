# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2021 INRIA.
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


from siconos.modeling import (
    NewtonImpactFrictionNSL,
    NonSmoothDynamicalSystem,
    interactions,
    copy_t,
)
from siconos.integrators import MoreauJeanOSI
from siconos.nonsmooth_formulations import FrictionContact
from siconos.simulation import TimeStepping, TimeDiscretisation
import siconos.numerics as sn

from siconos.mechanics.collision.bullet import SiconosBulletCollisionManager

from siconos.mechanics.collision import (
    SiconosBox,
    SiconosPlane,
    RigidBodyDS,
    SiconosContactor,
    SiconosContactorSet,
)

from numpy import zeros
from numpy.linalg import norm
import numpy as np
import siconos.plot_config as sicoplot

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)

t0 = 0  # start time
T = 20  # end time
h = 0.005  # time step

g = 9.81  # gravity

theta = 0.5  # theta scheme

#
# dynamical system
#
position_init = 10
velocity_init = 0

initial_position = np.array([0, 0, position_init, 1, 0, 0, 0], dtype=np.float64)
initial_velocity = np.array([0, 0, velocity_init, 0, 0, 0], dtype=np.float64)
inertia = np.eye(3, dtype=np.float64, order="F")
mass = 1.0
weight = np.array([0, 0, -mass * g], dtype=np.float64)


def makeBox(pos=initial_position, vel=initial_velocity):
    box = SiconosBox(1.0, 1.0, 1.0)

    # A Bullet Dynamical System : a shape + a mass (1.0) + position and velocity
    body = RigidBodyDS(initial_position, initial_velocity, mass, inertia)

    # set external forces
    body.setConstantFext(weight, copy_t)

    # Add the shape, wrapped in a SiconosContactor, to the body's
    # contactor set.
    body.contactors().append(SiconosContactor(box))

    return body


# Initial box body
body = makeBox()

# set external forces
weight = np.array([0, 0, -body.scalarMass * g], dtype=np.float64)
body.setConstantFext(weight, copy_t)

#
# Model
#
bouncingBox = NonSmoothDynamicalSystem(t0, T)

# add the dynamical system to the non smooth dynamical system
bouncingBox.insertDynamicalSystem(body)

#
# Simulation
#

# (1) OneStepIntegrators
osi = MoreauJeanOSI(theta)

ground = SiconosPlane()
groundOffset = np.array([0, 0, -0.5, 1, 0, 0, 0], dtype=np.float64)

# (2) Time discretisation --
timedisc = TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = FrictionContact(3)

osnspb.numericsSolverOptions().iparam[0] = 1000
osnspb.numericsSolverOptions().dparam[0] = 1e-5
osnspb.setMaxSize(16384)
osnspb.setMStorageType(sn.params.NM_SPARSE_BLOCK)
osnspb.setNumericsVerboseMode(False)

# keep previous solution
osnspb.setKeepLambdaAndYState(True)


# (4) non smooth law
nslaw = NewtonImpactFrictionNSL(0.8, 0.0, 0.0, 3)

# (5) broadphase contact detection
broadphase = SiconosBulletCollisionManager()

# insert a non smooth law for contactors id 0
broadphase.insertNonSmoothLaw(nslaw, 0, 0)

# The ground is a static object
# we give it a group contactor id : 0
scs = SiconosContactorSet()
scs.append(SiconosContactor(ground))
broadphase.addStaticBody(scs, groundOffset)

# (6) Simulation setup with (1) (2) (3) (4) (5)
simulation = TimeStepping(bouncingBox, timedisc)
simulation.insertInteractionManager(broadphase)

simulation.insertIntegrator(osi)
simulation.insertNonSmoothProblem(osnspb)


# Get the values to be plotted
# ->saved in a matrix dataPlot

N = int((T - t0) / h)
dataPlot = zeros((N + 1, 4))

#
# numpy pointers on dense Siconos vectors
#
q = body.q()
v = body.twist()

#
# initial data
#
dataPlot[0, 0] = t0
dataPlot[0, 1] = q[2]
dataPlot[0, 2] = v[2]

k = 1

# time loop
new_position = np.array([0, 0, 3, 1, 0, 0, 0], dtype=np.float64)
new_velocity = np.array([0, 0, 1, 0, 0, 0], dtype=np.float64)

while simulation.hasNextEvent():

    # Add a second box dynamically to the simulation
    if k == 100:
        ds = makeBox(new_position, new_velocity)
        bouncingBox.insertDynamicalSystem(ds)

    simulation.computeOneStep()

    dataPlot[k, 0] = simulation.nextTime()
    dataPlot[k, 1] = q[2]
    dataPlot[k, 2] = v[2]

    # if (broadphase.collisionWorld().getDispatcher().getNumManifolds() > 0):
    if (
        broadphase.statistics().new_interactions_created
        + broadphase.statistics().existing_interactions_processed
    ) > 0:
        if bouncingBox.topology().numberOfIndexSet() == 2:
            index1 = interactions(simulation.indexSet(1))
            if len(index1) == 4:
                dataPlot[k, 3] = (
                    norm(index1[0].lambda_python(1))
                    + norm(index1[1].lambda_python(1))
                    + norm(index1[2].lambda_python(1))
                    + norm(index1[3].lambda_python(1))
                )

    k += 1
    simulation.nextStep()

#
# comparison with the reference file
#
ref = np.loadtxt("result_dynamic.ref", skiprows=1)

print("norm(dataPlot - ref) = {0}".format(norm(dataPlot - ref)))


if norm(dataPlot - ref) > 1e-11:
    print("Warning. The result is rather different from the reference file.")


#
# plots
#

if enable_plot:
    plt.subplot(511)
    plt.title("position")
    plt.plot(dataPlot[0:k, 0], dataPlot[0:k, 1])
    y = plt.ylim()
    plt.plot(ref[0:k, 0], ref[0:k, 1])
    plt.ylim(y)
    plt.grid()
    plt.subplot(513)
    plt.title("velocity")
    plt.plot(dataPlot[0:k, 0], dataPlot[0:k, 2])
    y = plt.ylim()
    plt.plot(ref[0:k, 0], ref[0:k, 2])
    plt.ylim(y)
    plt.grid()
    plt.subplot(515)
    plt.plot(dataPlot[0:k, 0], dataPlot[0:k, 3])
    y = plt.ylim()
    plt.plot(ref[0:k, 0], ref[0:k, 3])
    plt.ylim(y)
    plt.title("lambda")
    plt.grid()
    plt.savefig("result_dynamic.png")
    if enable_plot:
        plt.show()

np.savetxt("BouncingBoxDynamic-py.dat", dataPlot)
