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

/*!\file PrismaticTest.cpp
  \brief \ref

  Simulation with a Time-Stepping scheme.
*/

#include <KneeJointR.hpp>
#include <PrismaticJointR.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>

#include "GeomTools.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;
    unsigned int qDim = 7;
    unsigned int nDim = 6;
    double t0 = 0;     // initial computation time
    double h = 0.001;  // time step
    double T = 10;
    double theta = 1.0;  // theta for MoreauJeanOSI integrator
    double g = 9.81;     // Gravity
    double m = 1.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ...\n";

    // -- Initial positions and velocities --
    auto q10 = std::make_shared<Vector>(qDim);
    auto v10 = std::make_shared<Vector>(nDim);
    auto I1 = std::make_shared<Matrix>(3, 3);
    v10->zero();
    I1->eye();
    q10->zero();
    q10->setValue(0, 1);
    q10->setValue(1, 1);
    q10->setValue(2, 1);

    double angle = numbers::pi / 5;
    Vector V1(3);
    V1.zero();
    V1.setValue(0, 3);
    V1.setValue(1, 2);
    V1.setValue(2, 1);
    double Vnorm = V1.norm2();
    V1.setValue(0, V1.getValue(0) / Vnorm);
    V1.setValue(1, V1.getValue(1) / Vnorm);
    V1.setValue(2, V1.getValue(2) / Vnorm);
    q10->setValue(3, cos(angle));
    q10->setValue(4, V1.getValue(0) * sin(angle));
    q10->setValue(5, V1.getValue(1) * sin(angle));
    q10->setValue(6, V1.getValue(2) * sin(angle));

    // -- The dynamical system --
    auto beam1 = std::make_shared<siconos::modeling::NewtonEulerDS>(q10, v10, m, I1);
    // -- Set external forces (weight) --
    Vector weight{nDof};
    weight.setZero();
    weight(2) = -m * g;
    beam1->etConstantFExt(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // Interaction ball-floor
    // -- prismatic axis 0,0,1 in absolute frame: ball can only move in Z
    auto axis1 = std::make_shared<Vector>(3);
    axis1->setValue(0, 0);
    axis1->setValue(1, 0);
    axis1->setValue(2, 1);

    auto relation1 = std::make_shared<siconos::joints::PrismaticJointR>(axis1, true, beam1);

    auto H1 = std::make_shared<Matrix>(relation1->numberOfConstraints(), qDim);
    H1->zero();
    relation1->setJachq(H1);

    auto nslaw1 = std::make_shared<siconos::modeling::EqualityConditionNSL>(
        relation1->numberOfConstraints());

    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);
    // -------------
    // --- Model ---
    // -------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    bouncingBall->insertDynamicalSystem(beam1);
    bouncingBall->link(inter1, beam1);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- OneStepIntegrators --
    auto OSI1 = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

    // -- OneStepNsProblem --
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::Equality>();

    auto s =
        std::make_shared<siconos::simulation::TimeStepping>(bouncingBall, t, OSI1, osnspb);
    s->setNewtonTolerance(5e-4);
    s->setNewtonMaxIteration(50);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    int N = 2000;  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 8;
    Matrix dataPlot(N, outputSize);

    auto q1 = beam1->q();
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 0;

    auto start = std::chrono::system_clock::now();
    std::vector<unsigned int> dimIndex(2);
    decltype(dimIndex) startIndex(4);
    int cmp = 0;
    for (cmp = 0; cmp < N; cmp++) {
      // solve ...
      s->advanceToEvent();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q1)(0);
      dataPlot(k, 2) = (*q1)(1);
      dataPlot(k, 3) = (*q1)(2);
      dataPlot(k, 4) = (*q1)(3);
      dataPlot(k, 5) = (*q1)(4);
      dataPlot(k, 6) = (*q1)(5);
      dataPlot(k, 7) = (*q1)(6);

      s->nextStep();

      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("PrismaticTest.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-11;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "PrismaticTest.ref", eps)) >
        eps)
      return 1;

    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
