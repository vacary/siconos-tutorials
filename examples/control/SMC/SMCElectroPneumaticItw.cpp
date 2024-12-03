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

/* !\file SMCElectroPneumatic.cpp
  \brief Simulation of an Electropneumatic setup controlled with a Twisting algorithm
  O. Huber
  */

#include <SiconosControl.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <string>
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;
using namespace std;

// main program
int main(int argc, char* argv[]) {
  // User-defined parameters
  unsigned ndof = 4;         // Number of degrees of freedom of your system
  double t0 = 0.0;           // Starting time
  double T = 100.0;          // Total simulation time
  double h = 1.0e-3;         // Time step for simulation
  double hControl = 1.0e-1;  // Time step for control
  double y0 = 0.0;
  double v0 = 0.0;
  double pP0 = 501570.871800124;
  double pN0 = 501570.871800124;

  auto param = std::make_shared<Vector>(3);
  (*param)(0) = 1e-2;
  (*param)(1) = 2. / 3.;
  (*param)(2) = 10.;
  // ================= Creation of the model =======================
  // Steps:
  // - create a Dynamical System
  // - add a Simulation to the model

  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  // First System:

  // Matrix declaration
  // For the DynamycalSystem
  auto x0 = std::make_shared<Vector>(ndof);
  (*x0)(0) = pP0;
  (*x0)(1) = pN0;
  (*x0)(2) = v0;
  (*x0)(3) = y0;

  // Dynamical Systems
  auto plant= std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(x0, "electro_pneumatic_nantesPlugin:computef", "electro_pneumatic_nantesPlugin:computeJacfx");
  plant->setzPtr(param);

  // -------------
  // --- Model process ---
  // -------------
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(t0, T, h);

  simLsodar->addDynamicalSystem(plant, "plant");
  // ------------------
  // --- Simulation ---
  // ------------------
  // TimeDiscretisation
  // Control stuff
  // For the Sensor
  auto sensorC = std::make_shared<Matrix>(4, 4);
  sensorC->setIdentity();
  auto sensor = std::make_shared<siconos::control::LinearSensor>(plant, sensorC);
  // add the sliding mode controller
  auto twisting = std::make_shared<siconos::control::Twisting>(sensor, hControl);

  twisting->seth("electro_pneumatic_nantesPlugin:computeh");
  twisting->setg("electro_pneumatic_nantesPlugin:computeg");
  twisting->setJachx("electro_pneumatic_nantesPlugin:computeJachx");
  twisting->setJachlambda("electro_pneumatic_nantesPlugin:computeJachlambda");
  twisting->setJacgx("electro_pneumatic_nantesPlugin:computeJacgx");
  twisting->setJacglambda("electro_pneumatic_nantesPlugin:computeJacglambda");

  simLsodar->addSensor(sensor, hControl);
  simLsodar->addActuator(twisting, hControl);

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Simulation initialisation ..." << endl << endl;
  // initialise the process and the ControlManager
  simLsodar->initialize();
  //  (std::static_pointer_cast<LsodarOSI>(simLsodar->integrator()))->setJT(1);
  //  (std::static_pointer_cast<LsodarOSI>(simLsodar->integrator()))->setMaxOrder(0, 5);
  simLsodar->run();
  auto& data = *simLsodar->data();
  siconos::algebra::io::write("SMCElectroPneumaticItw.dat", data,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);

  double error = 0.0, eps = 1e-12;
  if ((error = siconos::algebra::io::compareRefFile(data, "SMCElectroPneumaticItw.ref", eps)) >
      eps)
    return 1;
  else
    return 0;

  /*
    // Comparison with a reference file
    SiconosMatrix dataPlotRef(dataPlot);
    dataPlotRef.setZero();
    ioMatrix::read("SMCExampleImplicit.ref", "ascii", dataPlotRef);
    std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;
    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout << "Warning. The results is rather different from the reference file." <<
  std::endl; std::cout << (dataPlot - dataPlotRef).normInf() << std::endl; return 1;
    }
  #else
    return 0;
  */
  return 0;
}
