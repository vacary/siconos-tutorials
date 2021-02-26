#!/usr/bin/env python

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

import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import subplot, title, plot, grid, savefig
import numpy as np
import siconos.kernel as sk
import scipy.linalg as la

t0 = 0.0   # start time
T = 1.0     # end time
h = 1.0e-3   # time step
Vinit = 1.0
theta = 0.5
N = int((T-t0) / h)

# Dynamical systems
A = np.zeros((2, 2))
x0 = np.array([Vinit, -Vinit])
B = 2.0 * np.eye(2)
C = np.eye(2)
D = np.zeros((2,2))

process = sk.FirstOrderLinearDS(x0, A)

myNslaw = sk.RelayNSL(2)
myProcessRelation = sk.FirstOrderLinearR(C, B)
#myProcessRelation.setDPtr(D)

myProcessInteraction = sk.Interaction(myNslaw,
                                   myProcessRelation)

# NSDS
simplerelay = sk.NonSmoothDynamicalSystem(t0,T)
simplerelay.insertDynamicalSystem(process)
simplerelay.link(myProcessInteraction,process)

#myProcessRelation.computeJachx(0, x0, x0 , x0, C)

# Simulation
s = sk.TimeStepping(simplerelay,
                    sk.TimeDiscretisation(t0, h),
                    sk.EulerMoreauOSI(theta),
                    sk.Relay())
#s.setComputeResiduY(True)
#s.setComputeResiduR(True)

# matrix to save data
dataPlot = np.empty((N+1,7))
k=0
dataPlot[k, 0] = t0
dataPlot[k, 1:3] = process.x()
dataPlot[k, 3] = myProcessInteraction.lambda_(0)[0]
dataPlot[k, 4] = myProcessInteraction.lambda_(0)[1]
dataPlot[k, 5] = myProcessInteraction.y(0)[0]
dataPlot[k, 6] = myProcessInteraction.y(0)[1]
# time loop
k = 1
print('start computation')
while(s.hasNextEvent()):
     if not (k%50):
         sys.stdout.write('.')

     s.computeOneStep()
     dataPlot[k, 0] = s.nextTime()
     dataPlot[k, 1] = process.x()[0]
     dataPlot[k, 2] = process.x()[1]
     dataPlot[k, 3] = myProcessInteraction.lambda_(0)[0]
     dataPlot[k, 4] = myProcessInteraction.lambda_(0)[1]
     dataPlot[k, 5] = myProcessInteraction.y(0)[0]
     dataPlot[k, 6] = myProcessInteraction.y(0)[1]
     k += 1
     s.nextStep()
sys.stdout.write('\n')
# save to disk
np.savetxt('SimpleExampleRelay_py.dat', dataPlot)

dataRef = np.loadtxt('SimpleExampleRelay_py.ref')

print('Comparison with reference file -  error = ',la.norm(dataPlot-dataRef) )



assert(la.norm(dataPlot-dataRef) <= 1e-12)


# plot interesting stuff
subplot(311)
title('x_1')
plot(dataPlot[:,0], dataPlot[:,1])
grid()
subplot(312)
title('x_2')
plot(dataPlot[:,0], dataPlot[:,2])
grid()
subplot(313)
plot(dataPlot[:,0], dataPlot[:,3])
title('lambda')
grid()
savefig("SimpleRelay_py.png")
