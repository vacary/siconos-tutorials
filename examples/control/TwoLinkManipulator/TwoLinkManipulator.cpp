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

// =============================== Robot arm sample (HuMAnsPa10)
// ===============================
//
// see modelRobot1.jpg for complete system view.
//
// Keywords: LagrangianDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// =============================================================================================

#define _USE_MATH_DEFINES
#include <math.h>

#include <SiconosKernel.hpp>
#include <SiconosPointers.hpp>
#include <chrono>

#define PI 3.14159265

using namespace std;

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 2;  // degrees of freedom for robot arm
    double t0 = 0;          // initial computation time
    double T = 3.0;         // final computation time
    double h = 1e-4;        // time step
    double criterion = 1e-8;
    unsigned int maxIter = 20000;
    double e = 0.7;  // nslaw
    double e2 = 0.0;
    double L = 0.0;
    int test = 0;
    int nimpact = 0;

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // --- DS: manipulator arm ---

    // The dof are angles between ground and arm and between differents parts of the arm. (See
    // corresponding .pdf for more details)

    // Initial position (angles in radian)
    Vector q0(nDof), v0(nDof);
    q0.setZero();
    v0.setZero();
    q0(0) = 0.9;
    q0(1) = -1.6;
    auto z = std::make_shared<Vector>(nDof * 12);
    (*z)(0) = q0(0);
    (*z)(1) = q0(1);
    (*z)(2) = v0(0);
    (*z)(3) = v0(1);
    (*z)(4) = 0;
    (*z)(5) = 0;
    (*z)(6) = 0;
    (*z)(7) = 0;
    (*z)(8) = 0;
    (*z)(9) = 0;
    (*z)(10) = 0;
    (*z)(11) = PI;
    (*z)(12) = 0;
    (*z)(13) = 0;
    (*z)(14) = 0;
    (*z)(15) = 0;
    (*z)(16) = 0;
    (*z)(17) = 0;
    (*z)(22) = 0;
    (*z)(23) = 0;

    auto arm = std::make_shared<siconos::modeling::LagrangianDS>(
        siconos::pointers::createSPtr(q0), siconos::pointers::createSPtr(v0));

    // external plug-in
    arm->setComputeMassFunction("Two-linkPlugin", "mass");
    arm->setComputeFGyrFunction("Two-linkPlugin", "FGyr");
    arm->setComputeJacobianFGyrqDotFunction("Two-linkPlugin", "jacobianVFGyr");
    arm->setComputeJacobianFGyrqFunction("Two-linkPlugin", "jacobianFGyrq");
    arm->setComputeFIntFunction("Two-linkPlugin", "U");
    arm->setComputeJacobianFIntqDotFunction("Two-linkPlugin", "jacobFintV");
    arm->setComputeJacobianFIntqFunction("Two-linkPlugin", "jacobFintQ");
    arm->setzPtr(z);

    // -------------------
    // --- Interactions---
    // -------------------

    //  - one with Lagrangian non linear relation to define contact with ground
    //  Both with newton impact nslaw.
    // -- relations --

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    // auto relation=
    // std::make_shared<siconos::modeling::LagrangianScleronomousR>("Two-linkPlugin:h0",
    // "Two-linkPlugin:G0"); auto inter=
    // std::make_shared<siconos::modeling::Interaction>(nslaw, relation);
    auto relation01 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "Two-linkPlugin:h01", "Two-linkPlugin:G01");
    auto inter01 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation01);
    auto relation02 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "Two-linkPlugin:h02", "Two-linkPlugin:G02");
    auto inter02 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation02);

    auto H10 = std::make_shared<Matrix>(1, 2);
    auto b10 = std::make_shared<Vector>(1);
    H10->setZero();
    (*H10)(0, 0) = -1;
    (*b10)(0) = PI;

    auto nslaw2 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e2);
    auto relation10 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H10, *b10);
    auto inter10 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation10);

    auto H11 = std::make_shared<Matrix>(1, 2);
    auto b11 = std::make_shared<Vector>(1);
    H11->setZero();
    (*H11)(0, 0) = 1;
    (*b11)(0) = 0;

    auto relation11 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H11, *b11);
    auto inter11 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation11);

    auto H20 = std::make_shared<Matrix>(1, 2);
    auto b20 = std::make_shared<Vector>(1);
    H20->setZero();
    (*H20)(0, 1) = -1;
    (*b20)(0) = 0.0001;

    auto relation20 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H20, *b20);
    auto inter20 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation20);

    auto H21 = std::make_shared<Matrix>(1, 2);
    auto b21 = std::make_shared<Vector>(1);
    H21->setZero();
    (*H21)(0, 1) = 1;
    (*b21)(0) = PI - 0.0001;

    auto relation21 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H21, *b21);
    auto inter21 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation21);

    // -------------
    // --- Model ---
    // -------------

    auto Manipulator = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // add the dynamical system in the non smooth dynamical system
    Manipulator->insertDynamicalSystem(arm);

    // link the interaction and the dynamical system
    Manipulator->link(inter01, arm);
    Manipulator->link(inter02, arm);
    // link the interaction and the dynamical system
    Manipulator->link(inter10, arm);
    Manipulator->link(inter11, arm);
    // link the interaction and the dynamical system
    Manipulator->link(inter20, arm);
    Manipulator->link(inter21, arm);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    auto s = std::make_shared<siconos::simulation::TimeStepping>(Manipulator, t);
    s->setNewtonTolerance(criterion);
    s->setNewtonMaxIteration(maxIter);
    // -- OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(0.500001);
    s->insertIntegrator(OSI);

    // -- OneStepNsProblem --
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    s->insertNonSmoothProblem(osnspb);
    // OneStepNSProblem  osnspb=
    // std::make_shared<siconos::nonsmooth_formulations::LCP>(s,"name","Lemke",200001,
    // 0.00001);
    cout << "=== End of model loading === \n";

    // =========================== End of model definition ===========================

    // ================================= Computation
    unsigned int k = 0;
    unsigned int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 14;
    Matrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    auto q = arm->q();
    auto v = arm->velocity();
    auto p = arm->p(1);

    // Initialization of the dicrete parameter z needs the followinf first computations
    arm->computeJacobianFIntq(t0);
    arm->computeFint(*v, *q, t0);
    arm->computeJacobianFIntqDot(t0);
    inter01->computeOutput(t0, 0);
    inter02->computeOutput(t0, 0);
    inter10->computeOutput(t0, 0);
    inter20->computeOutput(t0, 0);
    inter11->computeOutput(t0, 0);
    inter21->computeOutput(t0, 0);

    // EventsManager * eventsManager = s->eventsManager();

    dataPlot(k, 0) = Manipulator->t0();
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*q)(1);
    dataPlot(k, 3) = (*inter02->y(0))(0);
    dataPlot(k, 4) = (*v)(0);
    dataPlot(k, 5) = (*v)(1);
    dataPlot(k, 6) = (*inter01->y(0))(0) - 2;
    dataPlot(k, 7) = nimpact;  //(*inter->y(1))(1);
    dataPlot(k, 8) = (*z)(6);
    dataPlot(k, 9) = L;  //(*z)(4);
    dataPlot(k, 10) = test;
    dataPlot(k, 11) = (*p)(1);
    dataPlot(k, 12) = (*z)(22);
    dataPlot(k, 13) = (*z)(23);

    cout << "====> Start computation ... " << endl << endl;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    while (k < N) {
      (*z)(0) = (*q)(0);
      (*z)(1) = (*q)(1);
      (*z)(2) = (*v)(0);
      (*z)(3) = (*v)(1);
      (*z)(16) = (*z)(14);
      (*z)(17) = (*z)(15);
      (*z)(20) = (*z)(18);
      (*z)(21) = (*z)(19);

      // get current time step
      k++;

      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*inter02->y(0))(0);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*inter01->y(0))(0) - 2;
      dataPlot(k, 7) = nimpact;  //(*inter->y(1))(1);
      dataPlot(k, 8) = (*z)(6);
      if (test == 3)
        dataPlot(k, 9) = (*z)(4) / h;
      else
        dataPlot(k, 9) = (*z)(4);
      dataPlot(k, 10) = test;
      dataPlot(k, 12) = (*z)(22);
      dataPlot(k, 13) = (*z)(23);

      s->advanceToEvent();
      dataPlot(k, 11) = (*p)(1);
      (*z)(4) = (inter02->getLambda(1))(0);
      s->nextStep();

      //    controller during impacts accumulation phase before the first impact
      if ((dataPlot(k, 3) <= 0.01) && (test == 0) && (dataPlot(k, 6) < 0.6)) {
        (*z)(8) = dataPlot(k, 0);
        (*z)(5) = 0.65 + 0.1 * cos(2 * PI * ((*z)(8)) / (*z)(11));
        (*z)(7) = (*z)(9);
        arm->setComputeFIntFunction("Two-linkPlugin", "U10");
        test = 1;
      }

      //  controller during impacts accumulation phase after the first impact
      if ((dataPlot(k, 11) > 0) && (test == 1)) {
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFIntFunction("Two-linkPlugin", "U11");
        test = 2;
      }
      if ((dataPlot(k, 11) > 0) && (test == 2)) nimpact = nimpact + 1;

      // controller during constraint-motion phase.
      if ((dataPlot(k, 11) > 0) && (test == 2) &&
          (dataPlot(k, 7) - dataPlot(k - 1, 7) == 1))  // && (fabs((*inter->y(1))(1))<1e-8))
      {
        L = dataPlot(k, 0) - (*z)(8);
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFIntFunction("Two-linkPlugin", "U2");
        test = 3;
        nimpact = 0;
      }

      // change of control law with a particular design of the desired trajectory that
      // guarantee the take-off
      if ((trunc((dataPlot(k, 0) + h) / (*z)(11)) > trunc((dataPlot(k, 0)) / (*z)(11))) &&
          (test == 3)) {
        (*z)(10) = dataPlot(k, 0) + h;
        (*z)(8) = (*z)(12);
        arm->setComputeFIntFunction("Two-linkPlugin", "U3");
        test = 4;
        L = 0;
      }

      //  controller during free-motion phase
      if (((*z)(13) - 0.1 >= 0) && (test == 4)) {
        arm->setComputeFIntFunction("Two-linkPlugin", "U");
        test = 0;
        (*z)(13) = 0;
      }
    }
    cout << endl << "End of computation - Number of iterations done: " << k << endl;
    cout << "Computation Time :\n";
    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Computation time : " << elapsed << " ms\n";

    // --- Output files ---
    // --- Output files ---
    dataPlot.resize(k, outputSize);

    siconos::algebra::io::write("TwoLinkManipulator.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "TwoLinkManipulator.ref",
                                                      eps)) > eps)
      return 1;
    else
      return 0;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
