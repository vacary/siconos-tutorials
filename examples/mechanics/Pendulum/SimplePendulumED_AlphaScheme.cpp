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

// =============================== Double Pendulum Example ===============================
//
// Author: Vincent Acary
//
// Keywords: LagrangianDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// =============================================================================================

#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>

using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

constexpr double gravity = 9.8100;

// User-defined main parameters
unsigned int nDof = 2;                 // degrees of freedom for robot arm
double L = 1.0;                        // Length of the pendulum
double InitAngle = numbers::pi / 3.0;  // Initial inclination angle
double m = 1.0;                        // Mass of the pendulum
double t0 = 0;                         // initial computation time
double T = 10.0;                       // final computation time
double h = 0.001;                      // time step
unsigned int N = ceil(T / h) + 1;      // Number of points to be saved
double e = 0.9;                        // nslaw
double _rho = 0.99;
bool IsHandleVelConstraint = false;
int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;

    // --- DS: Simple Pendulum ---

    // Initial position (angles in radian)
    auto q0 = std::make_shared<Vector>(nDof);
    auto v0 = std::make_shared<Vector>(nDof);
    (*q0).zero();
    (*v0).zero();
    (*q0)(0) = L * sin(InitAngle);
    (*q0)(1) = L * cos(InitAngle);

    auto Mass = std::make_shared<Matrix>(nDof, nDof);
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    auto simplependulum =
        std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, v0, Mass);

    std::vector<double> zparams = {L};
    auto zz = std::make_shared<Vector>(zparams);
    simplependulum->setzPtr(zz);

    auto ForceExtern = std::make_shared<Vector>(nDof);
    (*ForceExtern)(0) = 0.0;
    (*ForceExtern)(1) = m * gravity;
    simplependulum->setFExtPtr(ForceExtern);

    // -------------------
    // --- Interactions---
    // -------------------

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SimplePendulumBilateralConstraintPlugin:h0",
        "SimplePendulumBilateralConstraintPlugin:G0",
        "SimplePendulumBilateralConstraintPlugin:G0dot");
    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // -------------
    // --- Model ---
    // -------------

    auto Pendulum = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    Pendulum->insertDynamicalSystem(simplependulum);
    Pendulum->link(inter, simplependulum);
    // ----------------
    // --- Simulation ---
    // ----------------

    // 1. Time discretization
    auto TimeDiscret = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // 2. Integration solver for one step
    auto OSI =
        std::make_shared<siconos::integrators::NewMarkAlphaOSI>(_rho, IsHandleVelConstraint);
    // 3. Nonsmooth problem
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto acceleration = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto position = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    // 4. Simulation with (1), (2), (3)
    auto EDscheme = std::make_shared<siconos::simulation::EventDriven>(Pendulum, TimeDiscret);
    EDscheme->insertIntegrator(OSI);
    EDscheme->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_ED_IMPACT);
    EDscheme->insertNonSmoothProblem(acceleration,
                                     siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC);
    EDscheme->insertNonSmoothProblem(position,
                                     siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_POS);
    EDscheme->setPrintStat(true);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    auto eventsManager =
        EDscheme->eventsManager();            // ponters point to the "eventsManager" object
    auto _q = simplependulum->q();            // pointer points to the position vector
    auto _qdot = simplependulum->velocity();  // pointer points to the velocity
    simplependulum->initRhs(t0);
    simplependulum->computeRhs(t0);
    auto _qddot = simplependulum->acceleration();
    auto _g = inter->y(0);
    auto indexSet0 = Pendulum->topology()->indexSet(0);
    std::cout << "Size of IndexSet0: " << indexSet0->size() << endl;
    //-------------------- Save the output during simulation
    //---------------------------------------------------------
    Matrix DataPlot(N, 10);
    //------------- At the initial time
    //-----------------------------------------------------------------------------
    DataPlot(0, 0) = Pendulum->t0();
    DataPlot(0, 1) = (*_q)(0);      // Position X
    DataPlot(0, 2) = (*_q)(1);      // Position Y
    DataPlot(0, 3) = (*_qdot)(0);   // Velocity Vx
    DataPlot(0, 4) = (*_qdot)(1);   // Velocity Vy
    DataPlot(0, 5) = (*_qddot)(0);  // Acceleration ax
    DataPlot(0, 6) = (*_qddot)(1);  // Acceleration ay
    DataPlot(0, 7) = (*_g)(0);      // Contraint in position
    DataPlot(0, 8) = 0.0;           // Constraint in velocity
    DataPlot(0, 9) = 0.0;           // Reaction force

    //----------------------------------- Simulation starts
    //----------------------------------------------------------
    std::cout << "====> Start computation ... " << endl << endl;
    bool NSEvent = false;
    unsigned int NumberNSEvent = 0;
    unsigned int k = 0;

    auto start = std::chrono::system_clock::now();
    while ((EDscheme->hasNextEvent()) && (k < N)) {
      // std::cout << "--> k = " << k << std::endl;
      EDscheme->advanceToEvent();  // lead the simulation run from one event to the next
      //---------- detect the statue of the current event ------------------------------------
      if (eventsManager->nextEvent()->getType() ==
          siconos::simulation::EventType::NS)  // the current event is non-smooth
      {
        NSEvent = true;
      };
      EDscheme->processEvents();  // process the current event
      //------------------- get data at the beginning of non-smooth events
      //---------------------------
      auto _gdot = inter->y(1);
      auto _lambda = inter->lambda(2);

      if (NSEvent) {
        DataPlot(k, 0) = EDscheme->startingTime();  // instant at non-smooth event
        const auto& _qMemory = simplependulum->qMemory().getSiconosVector(0);
        const auto& _qdotMemory = simplependulum->velocityMemory().getSiconosVector(0);
        DataPlot(k, 1) = _qMemory(0);
        DataPlot(k, 2) = _qMemory(1);
        DataPlot(k, 3) = _qdotMemory(0);
        DataPlot(k, 4) = _qdotMemory(1);
        k++;
        ++NumberNSEvent;

        NSEvent = false;  // The next event is maybe smooth
      };
      //-------------------- get data at smooth events or at the end of non-smooth events
      //---------------
      DataPlot(k, 0) = EDscheme->startingTime();
      DataPlot(k, 1) = (*_q)(0);       // Position X
      DataPlot(k, 2) = (*_q)(1);       // Position Y
      DataPlot(k, 3) = (*_qdot)(0);    // Velocity Vx
      DataPlot(k, 4) = (*_qdot)(1);    // Velocity Vy
      DataPlot(k, 5) = (*_qddot)(0);   // Acceleration ax
      DataPlot(k, 6) = (*_qddot)(1);   // Acceleration ay
      DataPlot(k, 7) = (*_g)(0);       // Contraint in position
      DataPlot(k, 8) = (*_gdot)(0);    // Constraint in velocity
      DataPlot(k, 9) = (*_lambda)(0);  // Reaction force
      // go to the next time step
      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("SimplePendulumED_AlphaScheme.dat", DataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(
             DataPlot, "SimplePendulumED_AlphaScheme.ref", eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
