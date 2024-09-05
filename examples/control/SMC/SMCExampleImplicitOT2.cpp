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

/* !\file SMCExampleImplicitOT2.cpp
  \brief Two independent systems of dimension one controlled to slide
  on \f$x = 0\f$. The scheme from Su, Drakunov and Özgüner (2000) is used
  O. Huber
  */

#include <SiconosControl.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <string>
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;
using namespace std;

class MyDS : public siconos::modeling::FirstOrderLinearDS {
 public:
  MyDS(auto x0, auto A) : FirstOrderLinearDS(x0, A) {
    _b = std::make_shared<Vector>(x0->size());
  };
  void computeb(double time) override {
    double t = sin(50 * time);
    _b->setValue(0, t);
    _b->setValue(1, -t);
  };
};
// main program
int main(int argc, char* argv[]) {
  // User-defined parameters
  unsigned int ndof = 2;     // Number of degrees of freedom of your system
  double t0 = 0.0;           // Starting time
  double T = 1;              // Total simulation time
  double h = 1.0e-4;         // Time step for simulation
  double hControl = 1.0e-2;  // Time step for control
  double Xinit = 1.0;

  if (h > hControl) {
    THROW_EXCEPTION("hControl must be bigger than h");
  }

  // ================= Creation of the model =======================
  // Steps:
  // - create a Dynamical System
  // - add a Simulation to the model

  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  // First System:
  // dx/dt = Ax + u(t) + r
  // x(0) = x0
  // Note: r = Blambda, B defines in relation below.

  // Matrix declaration
  auto A = std::make_shared<Matrix>(ndof, ndof, 0);
  auto x0 = std::make_shared<Vector>(ndof);
  (*x0)(0) = Xinit;
  (*x0)(1) = -Xinit;
  auto sensorC = std::make_shared<Matrix>(2, 2);
  sensorC->eye();
  auto sensorD = std::make_shared<Matrix>(2, 2, 0);
  auto Csurface = std::make_shared<Matrix>(1, 2, 0);
  (*Csurface)(0, 1) = 1;
  auto Brel = std::make_shared<Matrix>(2, 1, 0);
  (*Brel)(1, 0) = 2;

  // Dynamical Systems
  auto processDS = std::make_shared<MyDS>(x0, A);
  // -------------
  // --- Model process ---
  // -------------
  auto sim = std::make_shared<siconos::control::ControlZOHSimulation>(t0, T, h);
  // sim->setSaveOnlyMainSimulation(true);
  sim->addDynamicalSystem(processDS);

  // ------------------
  // --- Simulation ---
  // ------------------
  // Control stuff
  // use a controlSensor
  auto sens = std::make_shared<siconos::control::LinearSensor>(processDS, sensorC);
  sim->addSensor(sens, hControl);
  // add the sliding mode controller
  auto act = std::make_shared<siconos::control::LinearSMCOT2>(sens);
  act->setCsurface(Csurface);
  act->setB(Brel);
  sim->addActuator(act, hControl);
  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Simulation initialisation ..." << endl << endl;
  // initialise the process and the ControlManager
  sim->initialize();

  // ==== Simulation loop =====
  cout << "====> Start computation ... " << endl << endl;
  sim->run();
  // --- Output files ---
  cout << "====> Output file writing ...\n";
  auto& dataPlot = *sim->data();

  siconos::algebra::io::write("SMCExampleImplicitOT2.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);

  double error = 0.0, eps = 1e-12;
  if ((error = siconos::algebra::io::compareRefFile(dataPlot, "SMCExampleImplicitOT2.ref",
                                                    eps)) > eps)
    return 1;
  else
    return 0;
}
