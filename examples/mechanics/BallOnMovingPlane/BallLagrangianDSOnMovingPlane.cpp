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

/*!\file BallOnMovingPlane.cpp
  \brief \ref EMBallOnMovingPlane - C++ input file, Time-Stepping version -
  V. Acary

  A Ball bouncing on a moving plane.
  This example shows how some precribed boundary conditions can be imposed.
  Simulation with a Time-Stepping scheme.
*/

#include <SiconosKernel.hpp>
#include <chrono>

using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;
using namespace std;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;       // degrees of freedom for the ball
    double t0 = 0;               // initial computation time
    double T = 10.0;             // final computation time
    double h = 0.005;            // time step
    double position_init = 1.0;  // initial position for lowest bead.
    double velocity_init = 0.0;  // initial velocity for lowest bead.
    double theta = 0.5;          // theta for MoreauJeanOSI integrator
    double R = 0.1;              // Ball radius
    double m = 1;                // Ball mass
    double g = 9.81;             // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    std::cout << "====> Model loading ...\n";

    auto Mass = std::make_shared<Matrix>(nDof, nDof);
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 2. / 5 * m * R * R;

    // -- Initial positions and velocities --
    auto q0 = std::make_shared<Vector>(nDof);
    auto v0 = std::make_shared<Vector>(nDof);
    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;

    // -- The dynamical system --
    auto ball = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, v0, Mass);

    // -- Set external forces (weight) --
    auto weight = std::make_shared<Vector>(nDof);
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

    // -- Moving Plane --

    // -- Initial positions and velocities --
    auto q02 = std::make_shared<Vector>(nDof);
    auto v02 = std::make_shared<Vector>(nDof);
    (*q02)(0) = 0.0;
    (*v02)(0) = -velocity_init;

    // -- The dynamical system --
    auto movingplane = std::make_shared<siconos::modeling::LagrangianDS>(q02, v02, Mass);

    // -- Set external forces (weight) --
    movingplane->setFExtPtr(weight);

    auto bd = std::make_shared<siconos::modeling::BoundaryCondition>(siconos::modeling::BoundaryCondition::Indices{0});
    bd->setComputePrescribedVelocityFunction("BallOnMovingPlanePlugin", "prescribedvelocity");

    movingplane->setBoundaryConditions(bd);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-plane
    //
    auto H = std::make_shared<Matrix>(1, 2 * nDof);
    (*H)(0, 0) = 1.0;
    (*H)(0, 3) = -1.0;

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H);

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // -------------
    // --- Model ---
    // -------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);
    bouncingBall->insertDynamicalSystem(movingplane);

    // link the interaction and the dynamical systems
    bouncingBall->link(inter, ball, movingplane);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(bouncingBall, t, OSI, osnspb);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 12;
    Matrix dataPlot(N + 1, outputSize);

    auto q = ball->q();
    auto v = ball->velocity();
    auto p = ball->p(1);
    auto qplane = movingplane->q();
    auto vplane = movingplane->velocity();
    auto pplane = movingplane->p(1);
    auto lambda = inter->lambda(1);
    auto y = inter->y(0);

    auto reaction = movingplane->reactionToBoundaryConditions();

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = (*q)(0);
    dataPlot(0, 6) = (*v)(0);
    dataPlot(0, 7) = (*qplane)(0);
    dataPlot(0, 8) = (*vplane)(0);
    dataPlot(0, 9) = (*qplane)(1);
    dataPlot(0, 10) = (*vplane)(1);
    dataPlot(0, 11) = (*reaction)(0);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();

    while (s->hasNextEvent() && k < 5000) {
      s->computeOneStep();
      // std::cout << "y :"<< std::endl;
      // y->display();
      // osnspb->display();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q)(1);
      dataPlot(k, 6) = (*v)(1);
      dataPlot(k, 7) = (*qplane)(0);
      dataPlot(k, 8) = (*vplane)(0);
      dataPlot(k, 9) = (*qplane)(1);
      dataPlot(k, 10) = (*vplane)(1);
      dataPlot(k, 11) = (*reaction)(0);

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
    siconos::algebra::io::write("result.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // Comparison with a reference file
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "BallOnMovingPlane.ref",
                                                      eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
