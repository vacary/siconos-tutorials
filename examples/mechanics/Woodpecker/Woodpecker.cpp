/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2023 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <SiconosKernel.hpp>
#include <chrono>

#include "WoodPeckerConsts.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  try {
    // ================= Model definition =================

    // User-defined main parameters
    unsigned int nDof = 3;  // degrees of freedom
    double t0 = 0;          // initial computation time
    double T = 0.3;         // final computation time
    double h = 0.00002;     // time step
    double theta = 0.5;     // theta for MoreauJeanOSI integrator;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    Matrix mass{nDof, nDof};
    mass(0, 0) = m_S + m_M;
    mass(0, 1) = m_S * l_M;
    mass(0, 2) = m_S * l_G;
    mass(1, 0) = m_S * l_M;
    mass(1, 1) = J_M + m_S * l_M * l_M;
    mass(1, 2) = m_S * l_M * l_G;
    mass(2, 0) = m_S * l_G;
    mass(2, 1) = m_S * l_M * l_G;
    mass(2, 2) = J_S + m_S * l_G * l_G;
    Matrix K{nDof, nDof};
    K(1, 1) = c_phi;
    K(1, 2) = -c_phi;
    K(2, 1) = -c_phi;
    K(2, 2) = c_phi;

    // -- Initial positions and velocities --
    Vector q0{nDof};
    q0(0) = y_0;
    q0(1) = phi_M_0;
    q0(2) = phi_S_0;

    Vector velocity0{nDof};
    velocity0(0) = v_0;
    velocity0(1) = omega_M_0;
    velocity0(2) = omega_S_0;

    auto dynamicalSystem =
        std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, velocity0, mass);
    dynamicalSystem->setStiffnessMatrix(K);

    dynamicalSystem->setComputeFextFunction(
        [](double time, Eigen::Ref<siconos::algebra::MapVectorType> fext) {
          fext(0) = -(m_S + m_M) * g;
          fext(1) = -m_S * l_M * g;
          fext(2) = -m_S * l_G * g;
        });

    // --------------------
    // --- Interactions ---
    // --------------------

    auto H1 = std::make_shared<Matrix>(2, nDof);
    (*H1)(0, 0) = 0;
    (*H1)(0, 1) = 0;
    (*H1)(0, 2) = -h_S;
    (*H1)(1, 0) = 1;
    (*H1)(1, 1) = l_M;
    (*H1)(1, 2) = l_G - l_S;
    auto b1 = std::make_shared<Vector>(2);
    (*b1)(0) = l_M + l_G - l_S - r_O;
    (*b1)(1) = 0;
    auto nslaw1 = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(
        eps_N_1, eps_T_123, mu_123, 2);

    auto relation1 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H1, *b1);

    auto H2 = std::make_shared<Matrix>(2, nDof);
    (*H2)(0, 0) = 0;
    (*H2)(0, 1) = h_M;
    (*H2)(0, 2) = 0;
    (*H2)(1, 0) = 1;
    (*H2)(1, 1) = r_M;
    (*H2)(1, 2) = 0;
    auto H3 = std::make_shared<Matrix>(2, nDof);
    (*H3)(0, 0) = 0;
    (*H3)(0, 1) = -h_M;
    (*H3)(0, 2) = 0;
    (*H3)(1, 0) = 1;
    (*H3)(1, 1) = r_M;
    (*H3)(1, 2) = 0;
    auto b2 = std::make_shared<Vector>(2);
    (*b2)(0) = r_M - r_O;
    (*b2)(1) = 0;
    auto b3 = std::make_shared<Vector>(2);
    (*b3)(0) = r_M - r_O;
    (*b3)(1) = 0;

    auto nslaw23 = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(
        eps_N_23, eps_T_123, mu_123, 2);

    auto relation2 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H2, *b2);
    auto relation3 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H3, *b3);

    auto I1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);

    auto I2 = std::make_shared<siconos::modeling::Interaction>(nslaw23, relation2);

    auto I3 = std::make_shared<siconos::modeling::Interaction>(nslaw23, relation3);

    // -------------
    // --- Model ---
    // -------------

    auto model = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    model->insertDynamicalSystem(dynamicalSystem);
    model->link(I1, dynamicalSystem);
    model->link(I2, dynamicalSystem);
    model->link(I3, dynamicalSystem);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- OneStepIntegrators --
    auto vOSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(2);

    auto s = std::make_shared<siconos::simulation::TimeStepping>(model, t, vOSI, osnspb);

    std::cout << "=== End of model loading === \n";

    // ================= Computation =================

    int k = 0;
    int N = ceil((T - t0) / h);

    // --- Get the values to be plotted ---
    unsigned int outputSize = 7;
    Matrix dataPlot(N + 1, outputSize);
    dataPlot(k, 0) = t0;
    for (int i = 0; i < (int)nDof; i++) {
      dataPlot(k, 2 * i + 1) = (dynamicalSystem->q_read())(i);
      dataPlot(k, 2 * i + 2) = (dynamicalSystem->velocity_read())(i);
    }

    // --- Time loop ---
    std::cout << "Start computation ... \n";
    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      // get current time step
      k++;

      // solve ...
      s->computeOneStep();

      // get values
      dataPlot(k, 0) = s->nextTime();
      for (int i = 0; i < (int)nDof; i++) {
        dataPlot(k, 2 * i + 1) = (dynamicalSystem->q_read())(i);
        dataPlot(k, 2 * i + 2) = (dynamicalSystem->velocity_read())(i);
      }

      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1 << std::endl;
    std::cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    std::cout << "====> Output file writing ..." << std::endl;
    siconos::algebra::io::write("Woodpecker.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "Woodpecker.ref", eps)) > eps)
      return 1;

    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
  return 0;
}
