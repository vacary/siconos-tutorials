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

using namespace std;
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;       // degrees of freedom for the ball
    double t0 = 0;               // initial computation time
    double T = 10;               // final computation time
    double h = 0.005;            // time step
    double position_init = 1.0;  // initial position for lowest bead.
    double velocity_init = 0.0;  // initial velocity for lowest bead.
    double theta = 0.5;          // theta for MoreauJeanOSI integrator
    double R = 0.5;              // Ball radius
    double m = 1;                // Ball mass
    double g = 9.81;             // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl;

    Matrix mass{nDof, nDof};
    mass.setZero();
    mass(0, 0) = m;
    mass(1, 1) = m;
    mass(2, 2) = 2. / 5 * m * R * R;

    // -- Initial positions and velocities --
    Vector q0{nDof};
    q0.setZero();
    q0(0) = position_init;
    Vector v0{nDof};
    v0.setZero();
    v0(0) = velocity_init;


    // -- The dynamical system --
    auto ball = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, v0, mass);

    auto q01 = std::make_shared<Vector>(nDof);
    auto v01 = std::make_shared<Vector>(nDof);
    (*q01)(0) = position_init + 2 * R + 0.1;
    (*v01)(0) = velocity_init;

    auto ball1 = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q01, v01, mass);

    // -- Set external forces (weight) --
    Vector weight{nDof};
    weight(0) = -m * g;
    ball->setConstantFext(weight);
    ball1->setConstantFext(weight);
    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);

    auto relation = std::make_shared<siconos::modeling::Lagrangian2d1DR>();

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    auto relation1 = std::make_shared<siconos::modeling::Lagrangian2d1DR>();

    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation1);

    // -------------
    // --- Model ---
    // -------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);
    bouncingBall->insertDynamicalSystem(ball1);

    // link the interaction and the dynamical system
    bouncingBall->link(inter, ball);
    bouncingBall->link(inter1, ball1, ball);

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

    int N = ceil((T - t0) / h) + 1;  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    Matrix dataPlot(N + 1, outputSize);

    auto q = ball->q();
    auto v = ball->velocity();
    auto p = ball->p(1);
    auto lambda = inter->lambda(1);
    auto q1 = ball1->q();
    auto v1 = ball->velocity();
    auto p1 = ball->p(1);
    auto lambda1 = inter1->lambda(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = (*q1)(0);
    dataPlot(0, 6) = (*v1)(0);
    dataPlot(0, 7) = (*p1)(0);
    dataPlot(0, 8) = (*lambda1)(0);
    // --- Time loop ---
    cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    auto start = std::chrono::system_clock::now();
    auto rpc = relation->pc1();
    auto nnc = relation->nc();

    auto rpc1 = relation1->pc1();
    auto rpc2 = relation1->pc2();
    auto nnc1 = relation1->nc();

    while (s->hasNextEvent()) {
      // a fake contact detection
      (*rpc)(0) = -R + (*q)(0);
      (*rpc)(1) = (*q)(1);
      (*nnc)(0) = 1.0;
      (*nnc)(1) = 0.0;

      (*rpc1)(0) = -R + (*q1)(0);
      (*rpc1)(1) = (*q1)(1);

      (*rpc2)(0) = R + (*q)(0);
      ;
      (*rpc2)(1) = (*q)(1);
      (*nnc1)(0) = 1.0;
      (*nnc1)(1) = 0.0;

      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q1)(0);
      dataPlot(k, 6) = (*v1)(0);
      dataPlot(k, 7) = (*p1)(0);
      dataPlot(k, 8) = (*lambda1)(0);
      s->nextStep();

      k++;
    }
    cout << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time : " << endl;
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Computation time : " << elapsed << " ms\n";

    // --- Output files ---
    cout << "====> Output file writing ...\n";
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("Ball2D_kernel_only.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "Ball2D_kernel_only.ref",
                                                      eps)) > eps)
      return 1;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
