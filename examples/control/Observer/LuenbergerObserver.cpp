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

/*!\file LuenbergerObserver.cpp
  \brief Academic example of a Luenberger Observer
  O. Huber.

  The controlled plant is a double integrator
  */

#include <LinearSensor.hpp>
#include <PID.hpp>
#include <SiconosControl.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <string>

using namespace std;
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  // ================= Creation of the model =======================

  // User-defined main parameters
  unsigned int nDof = 2;  // degrees of freedom for the system
  double t0 = 0;          // initial computation time
  double T = 100;         // final computation time
  double h = 0.05;        // time step
  double hControl = 2 * h;
  double position_init = 10;   // initial position for lowest bead.
  double velocity_init = 0.0;  // initial velocity for lowest bead.
  double xFinal = 0;           // final value
  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  cout << "====> Model loading ..." << endl << endl;

  auto A = std::make_shared<Matrix>(nDof, nDof);
  A->zero();
  (*A)(0, 1) = 1;
  (*A)(1, 0) = -1;

  auto B = std::make_shared<Matrix>(nDof, 1);
  (*B)(nDof - 1, 0) = 1.0;
  // -- Initial positions and velocities --
  auto x0 = std::make_shared<Vector>(nDof);
  (*x0)(0) = position_init;
  (*x0)(1) = velocity_init;

  // -- The dynamical system --
  auto doubleIntegrator = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(x0, A);

  // -------------
  // --- Model ---
  // -------------
  auto sim = std::make_shared<siconos::control::ControlZOHSimulation>(t0, T, h);

  sim->addDynamicalSystem(doubleIntegrator);

  // ------------------
  // --- Simulation ---
  // ------------------

  // use a controlSensor
  auto C = std::make_shared<Matrix>(1, 2, 0);
  (*C)(0, 0) = 1;
  auto sens = std::make_shared<siconos::control::LinearSensor>(doubleIntegrator, C);
  sim->addSensor(sens, hControl);

  // add the Observer
  auto L = std::make_shared<Matrix>(2, 1);
  (*L)(0, 0) = -8.71455866;
  (*L)(1, 0) = -24.02084375;
  auto xHat0 = std::make_shared<Vector>(2);
  (*xHat0)(0) = -10;
  (*xHat0)(1) = -5;
  auto obs = std::make_shared<siconos::control::LuenbergerObserver>(sens, *xHat0, C, L);
  sim->addObserver(obs, hControl);
  // add the PID controller
  auto K = std::make_shared<Vector>(3, 0);
  (*K)(0) = .25;
  (*K)(1) = .125;
  (*K)(2) = 2;
  auto act = std::make_shared<siconos::control::PID>(sens);
  act->setRef(xFinal);
  act->setK(K);
  act->setB(B);
  sim->addActuator(act, hControl);

  cout << "=== End of model loading === \n";
  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Initialisation ..." << endl << endl;
  // Initialize the model and the controlManager
  sim->initialize();
  // --- Time loop ---
  cout << "====> Start computation ... " << endl << endl;
  // ==== Simulation loop - Writing without explicit event handling =====
  sim->run();
  cout << endl << sim->dataLegend() << endl;

  // --- Output files ---
  cout << "====> Output file writing ...\n";
  auto& dataPlot = *sim->data();
  siconos::algebra::io::write("LuenbergerObserver.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);

  double error = 0.0, eps = 1e-12;
  if ((error = siconos::algebra::io::compareRefFile(dataPlot, "LuenbergerObserver.ref", eps)) >
      eps)
    return 1;
  else
    return 0;
}
