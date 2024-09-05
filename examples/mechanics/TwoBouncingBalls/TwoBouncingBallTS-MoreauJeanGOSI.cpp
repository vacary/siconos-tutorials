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

/*!\file BouncingBallTS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;  // degrees of freedom for the ball
    double t0 = 0;          // initial computation time
    double T = 3.0;
    double h = 0.001;            // time step
    double position_init = 0.1;  // initial position for lowest bead.
    double velocity_init = 0.0;  // initial velocity for lowest bead.
    double theta = 0.5;          // theta for MoreauJeanOSI integrator
    double R = 0.1;              // Ball radius
    double m1 = 1;               // Ball mass
    double m2 = 2.0;             // Ball mass
    double g = 9.81;             // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl;

    auto Mass = std::make_shared<Matrix>(nDof, nDof);
    (*Mass)(0, 0) = m1;
    (*Mass)(1, 1) = m1;
    (*Mass)(2, 2) = 2. / 5 * m1 * R * R;

    auto Mass2 = std::make_shared<Matrix>(nDof, nDof);
    (*Mass2)(0, 0) = m2;
    (*Mass2)(1, 1) = m2;
    (*Mass2)(2, 2) = 2. / 5 * m2 * R * R;

    // -- Initial positions and velocities --
    auto q0 = std::make_shared<Vector>(nDof);
    auto v0 = std::make_shared<Vector>(nDof);
    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;

    auto q0_2 = std::make_shared<Vector>(nDof);
    auto v0_2 = std::make_shared<Vector>(nDof);
    (*q0_2)(0) = position_init + 2 * R + 0.001;
    (*v0_2)(0) = velocity_init;

    // -- The dynamical system --
    auto ball1 = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, v0, Mass);
    auto ball2 = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0_2, v0_2, Mass2);

    // -- Set external forces (weight) --
    auto weight = std::make_shared<Vector>(nDof);
    (*weight)(0) = -m1 * g;
    ball1->setFExtPtr(weight);

    auto weight2 = std::make_shared<Vector>(nDof);
    (*weight2)(0) = -m2 * g;
    ball2->setFExtPtr(weight2);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.95;

    // Interaction ball-floor
    //
    auto H = std::make_shared<Matrix>(nDof, nDof);
    (*H)(0, 0) = 1.0;
    (*H)(1, 1) = 1.0;
    (*H)(2, 2) = 1.0;

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(e, e, 0.6, 3);
    auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H);

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // Interaction ball-ball
    //
    auto H_bb = std::make_shared<Matrix>(nDof, 2 * nDof);
    (*H_bb)(0, 0) = -1.0;
    (*H_bb)(1, 1) = -1.0;
    (*H_bb)(2, 2) = -1.0;
    (*H_bb)(0, 3) = 1.0;
    (*H_bb)(1, 4) = 1.0;
    (*H_bb)(2, 5) = 1.0;

    auto b_bb = std::make_shared<Vector>(3);
    (*b_bb)(0) = -2 * R;
    (*b_bb)(1) = 0.0;
    (*b_bb)(2) = 0.0;

    auto relation_bb = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H_bb, b_bb);

    auto inter_bb = std::make_shared<siconos::modeling::Interaction>(nslaw, relation_bb);

    // -------------
    // --- Model ---
    // -------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball1);
    bouncingBall->insertDynamicalSystem(ball2);

    // link the interaction and the dynamical system
    bouncingBall->link(inter, ball1);
    bouncingBall->link(inter_bb, ball1, ball2);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanGOSI>(theta);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::GlobalFrictionContact>(
        3, SICONOS_GLOBAL_FRICTION_3D_NSGS_WR);
    // auto osnspb=
    // std::make_shared<siconos::nonsmooth_formulations::GlobalFrictionContact>(3,SICONOS_GLOBAL_FRICTION_3D_ADMM);

    // auto osnspb=
    // std::make_shared<siconos::nonsmooth_formulations::GlobalFrictionContact>(3,SICONOS_GLOBAL_FRICTION_3D_NSN_AC);
    // auto osnspb=
    // std::make_shared<siconos::nonsmooth_formulations::GlobalFrictionContact>(3);
    assert(osnspb->numericsSolverOptions());
    SolverOptions* options = osnspb->numericsSolverOptions().get();
    // solver_options_print(options);
    options->dparam[SICONOS_DPARAM_TOL] = 1e-13;
    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(bouncingBall, t, OSI, osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 8;
    Matrix dataPlot(N + 1, outputSize);

    auto q1 = ball1->q();
    auto v1 = ball1->velocity();
    auto p1 = ball1->p(1);
    auto lambda = inter->lambda(1);
    auto q2 = ball2->q();
    auto v2 = ball2->velocity();
    auto p2 = ball2->p(1);
    // auto lambda = inter->lambda(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q1)(0);
    dataPlot(0, 2) = (*v1)(0);
    dataPlot(0, 3) = (*p1)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = (*q2)(0);
    dataPlot(0, 6) = (*v2)(0);
    dataPlot(0, 7) = (*p2)(0);
    // --- Time loop ---
    cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      osnspb->setNumericsVerboseMode(0);

      // std::cout << "############################  time step = " << k <<std::endl;
      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q1)(0);
      dataPlot(k, 2) = (*v1)(0);
      dataPlot(k, 3) = (*p1)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q2)(0);
      dataPlot(k, 6) = (*v2)(0);
      dataPlot(k, 7) = (*p2)(0);
      // osnspb->display();
      s->nextStep();

      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("TwoBouncingBallTS-MoreauJeanGOSI.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "TwoBouncingBallTS-MoreauJeanGOSI.ref", eps)) > eps)
      return 1;

    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
