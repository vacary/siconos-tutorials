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
//-----------------------------------------------------------------------
//
// Colpitts Oscillator
//-----------------------------------------------------------------------
#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <string>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  double t0 = 0.0;
  double T = 100.0;        // Total simulation time
  double h_step = 1.0e-3;  // Time step
  double L = 0.1;          // inductance
  double C1 = 2.0;         // capacitance
  double C2 = 0.8;         // capacitance
  double Rc = 10.0;        // resistance
  double Re = 20.0;        // resistance
  // double Rb = 0.5;    // resistance
  double alphaF = 0.99;
  double alphaR = 0.015;
  double VCC = 10;
  double VEE = 20;
  std::string Modeltitle = "Colpitts";

  try {
    // --- Dynamical system specification ---
    auto init_state = std::make_shared<Vector>(3);
    init_state->setZero();
    //    init_state->setValue(1,-1.0);
    auto LS_A = std::make_shared<Matrix>(3, 3);
    LS_A->setZero();

    LS_A->setValue(0, 0, -1.0 / (Rc * C1));
    LS_A->setValue(1, 0, -1.0 / (Rc * C2));
    LS_A->setValue(2, 0, -1.0 / L);

    LS_A->setValue(0, 1, 1.0 / C1 * (-1.0 / Rc));
    LS_A->setValue(1, 1, -1.0 / C2 * (1.0 / Rc + 1.0 / Re));
    LS_A->setValue(2, 1, -1.0 / L);

    LS_A->setValue(0, 2, 1.0 / (C1));
    LS_A->setValue(1, 2, 1.0 / (C2));
    LS_A->setValue(2, 2, 0.0);

    auto LS_b = std::make_shared<Vector>(3);

    (*LS_b)(0) = VCC / (Rc * C1);
    (*LS_b)(1) = 1.0 / C2 * (VCC / Rc - VEE / Re);
    (*LS_b)(2) = VCC / L;

    auto LSCollpitts =
        std::make_shared<siconos::modeling::FirstOrderLinearDS>(*init_state, *LS_A, *LS_b);

    // --- Interaction between linear system and non smooth system ---
    auto Int_C = std::make_shared<Matrix>(2, 3);
    Int_C->setZero();
    (*Int_C)(0, 0) = 1.0;
    (*Int_C)(1, 0) = 0.0;

    (*Int_C)(0, 1) = 1.0;
    (*Int_C)(1, 1) = 1.0;

    (*Int_C)(0, 2) = 0.0;
    (*Int_C)(1, 2) = 0.0;

    // auto Int_D= std::make_shared<Matrix>(2, 2,0.0);
    //  (*Int_D)(0, 0) = 1.0 / Rvalue;
    //  (*Int_D)(0, 1) = 1.0 / Rvalue;
    //  (*Int_D)(0, 2) = -1.0;
    //  (*Int_D)(1, 0) = 1.0 / Rvalue;
    //  (*Int_D)(1, 1) = 1.0 / Rvalue;
    //  (*Int_D)(1, 3) = -1.0;
    //  (*Int_D)(2, 0) = 1.0;
    //  (*Int_D)(3, 1) = 1.0;

    auto Int_B = std::make_shared<Matrix>(3, 2);
    Int_B->setZero();
    (*Int_B)(0, 0) = 1.0 / C1;
    (*Int_B)(1, 0) = (1.0 - alphaR) / C2;
    (*Int_B)(2, 0) = 0.0;
    (*Int_B)(0, 1) = -alphaF / C1;
    (*Int_B)(1, 1) = (1.0 - alphaF) / C2;
    (*Int_B)(2, 1) = 0.0;

    auto LTIRCollpitts =
        std::make_shared<siconos::modeling::FirstOrderLinearTIR>(*Int_C, *Int_B);
    // LTIRCollpitts->setConstantD(*Int_D);

    auto nslaw = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(2);

    auto InterCollpitts =
        std::make_shared<siconos::modeling::Interaction>(nslaw, LTIRCollpitts);

    // --- Model creation ---
    auto Collpitts = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    Collpitts->setTitle(Modeltitle);
    // add the dynamical system in the non smooth dynamical system
    Collpitts->insertDynamicalSystem(LSCollpitts);
    // link the interaction and the dynamical system
    Collpitts->link(InterCollpitts, LSCollpitts);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    double theta = 0.5;
    auto aOSI = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);

    // -- (2) Time discretisation --
    auto aTiDisc = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h_step);

    // -- (3) Non smooth problem

    auto aLCP = std::make_shared<siconos::nonsmooth_formulations::LCP>(SICONOS_LCP_LEMKE);
    aLCP->numericsSolverOptions()->dparam[0] = 1e-08;

    // auto aLCP= std::make_shared<siconos::nonsmooth_formulations::LCP>(SICONOS_LCP_ENUM);
    // aLCP->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS]=1;  //
    // Multiple solutions 0 or 1
    // aLCP->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_SEED]=4;  // choice of
    // seeds for multiple solutions
    // aLCP->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_USE_DGELS]=1;  // LS for
    // enum aLCP->setNumericsVerboseMode(1);

    // -- (4) Simulation setup with (1) (2) (3)
    auto aTS =
        std::make_shared<siconos::simulation::TimeStepping>(Collpitts, aTiDisc, aOSI, aLCP);

    int k = 0;
    double h = aTS->timeStep();
    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    Matrix dataPlot{N, 8};

    auto x = LSCollpitts->x();
    auto y = InterCollpitts->y(0);
    auto lambda = InterCollpitts->lambda(0);

    // For the initial time step:
    // time
    dataPlot(k, 0) = t0;

    dataPlot(k, 1) = (*x)(0);
    dataPlot(k, 2) = (*x)(1);
    dataPlot(k, 3) = (*x)(2);

    dataPlot(k, 4) = (*y)(0);
    dataPlot(k, 5) = (*y)(1);

    dataPlot(k, 6) = (*lambda)(0);
    dataPlot(k, 7) = (*lambda)(1);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    // --- Time loop  ---
    for (k = 1; k < N; ++k) {
      // solve ...
      aTS->computeOneStep();
      //  siconos::algebra::print(*aLCP);
      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = aTS->nextTime();

      dataPlot(k, 1) = (*x)(0);
      dataPlot(k, 2) = (*x)(1);
      dataPlot(k, 3) = (*x)(2);

      dataPlot(k, 4) = (*y)(0);
      dataPlot(k, 5) = (*y)(1);

      dataPlot(k, 6) = (*lambda)(0);
      dataPlot(k, 7) = (*lambda)(1);

      aTS->nextStep();
    }

    // --- elapsed time computing ---
    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";

    // Number of time iterations
    std::cout << "Number of iterations done: " << k << "\n";

    // dataPlot (ascii) output
    siconos::algebra::io::write("Colpitts.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // double error = 0.0, eps = 1e-12;
    // if ((error = siconos::algebra::io::compareRefFile(dataPlot, "Colpitts.ref", eps)) > eps) {
    //   if ((error = siconos::algebra::io::compareRefFile(dataPlot, "Colpitts-sol2.ref", eps)) >
    //       eps)
    //     return 1;
    // }
    return 0;
  }
  // --- Exceptions handling ---
  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
