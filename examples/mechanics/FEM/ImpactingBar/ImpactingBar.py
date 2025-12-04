import numpy as np
import siconos.modeling as sm
import siconos.integrators
import siconos.simulation
import siconos.nonsmooth_formulations
import matplotlib.pyplot as plt


# User-defined main parameters
nDof = 500  # degrees of freedom for the beam
t0 = 1e-8  # initial computation time
T = 0.0015  # final computation time
h = 2e-7  # time step
position_init = 0.00005  # initial position
velocity_init = -0.1  # initial velocity
epsilon = 0.5  # 1e-1
theta = 1 / 2.0 + epsilon  # theta for MoreauJeanOSI integrator
# theta = 1.0
E = 210e9  # young Modulus
S = 0.000314  # Beam Section 1 cm  for the diameter
# S=0.1
beam_length = 1.0  # length of the  beam
elem_length = beam_length / nDof  # length of an element
rho = 7800.0  # specific mass
# rho=1.0
g = 9.81  # Gravity
g = 0.0


# mass= SiconosMatrix(nDof,nDof,SPARSE,nDof)
# stiffness= SiconosMatrix(nDof,nDof,SPARSE,nDof)
mass = np.zeros((nDof, nDof), dtype=np.float64, order="F")
stiffness = np.zeros((nDof, nDof), dtype=np.float64, order="F")

stiffness[0, 0] = 1.0 * E * S / elem_length
stiffness[0, 1] = -1.0 * E * S / elem_length
mass[0, 0] = 1 / 3.0 * rho * S * elem_length
mass[0, 1] = 1 / 6.0 * rho * S * elem_length

for i in range(1, nDof - 1):
    stiffness[i, i] = 2.0 * E * S / elem_length
    stiffness[i, i - 1] = -1.0 * E * S / elem_length
    stiffness[i, i + 1] = -1.0 * E * S / elem_length
    mass[i, i] = 2 / 3.0 * rho * S * elem_length
    mass[i, i - 1] = 1 / 6.0 * rho * S * elem_length
    mass[i, i + 1] = 1 / 6.0 * rho * S * elem_length


stiffness[nDof - 1, nDof - 2] = -1.0 * E * S / elem_length
stiffness[nDof - 1, nDof - 1] = 1.0 * E * S / elem_length
mass[nDof - 1, nDof - 2] = 1 / 6.0 * rho * S * elem_length
mass[nDof - 1, nDof - 1] = 1 / 3.0 * rho * S * elem_length


q0 = np.full((nDof), position_init)
v0 = np.full((nDof), velocity_init)

bar = sm.LagrangianLinearTIDS(q0, v0, mass)
bar.setStiffnessMatrix(stiffness)
# bar.display()

weight = np.full((nDof), -g * rho * S / elem_length)
bar.setConstantFext(weight)

e = 0.0

H = np.zeros((1, nDof))
H[0, 0] = 1.0

nslaw = sm.NewtonImpactNSL(e)
relation = sm.LagrangianLinearTIR(H)
inter = sm.Interaction(nslaw, relation)

# -------------
# --- Model ---
# -------------
impactingBar = sm.NonSmoothDynamicalSystem(t0, T)

# add the dynamical system in the non smooth dynamical system
impactingBar.insertDynamicalSystem(bar)

# link the interaction and the dynamical system
impactingBar.link(inter, bar)


# ------------------
# --- Simulation ---
# ------------------

# -- (1) OneStepIntegrators --
OSI = siconos.integrators.MoreauJeanOSI(theta, 0.5)

# -- (2) Time discretisation --
t = siconos.simulation.TimeDiscretisation(t0, h)

# -- (3) one step non smooth problem
osnspb = siconos.nonsmooth_formulations.LCP()

s = siconos.simulation.TimeStepping(impactingBar, t, OSI, osnspb)

k = 0

N = int((T - t0) / h)
dataPlot = np.zeros((N + 1, 5))

q = bar.q()
v = bar.velocity()
p = bar.p(1)
lambda_ = inter.lambda_python(1)

# time loop
while s.hasNextEvent():
    s.computeOneStep()
    dataPlot[k, 0] = s.nextTime()
    # print("time=", dataPlot[k, 0])
    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = v[0]
    dataPlot[k, 3] = p[0] / h
    dataPlot[k, 4] = lambda_[0]

    k += 1
    s.nextStep()

dataPlot.resize(k, 5)

fig_size = [14, 14]
plt.rcParams["figure.figsize"] = fig_size

plt.subplot(411)
plt.title("position")
plt.plot(dataPlot[:, 0], dataPlot[:, 1])
plt.grid()
plt.subplot(412)
plt.title("velocity")
plt.plot(dataPlot[:, 0], dataPlot[:, 2])
plt.grid()
plt.subplot(413)
plt.plot(dataPlot[:, 0], dataPlot[:, 3])
plt.title("reaction")
plt.grid()
plt.subplot(414)
plt.plot(dataPlot[:, 0], dataPlot[:, 4])
plt.title("lambda")
plt.grid()


plt.show()
