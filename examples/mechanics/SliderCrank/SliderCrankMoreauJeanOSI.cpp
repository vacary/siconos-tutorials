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

/*!\file
  C++ input file, MoreauJeanOSI-Time-Stepping version
  T. Schindler, V. Acary

  Slider-crank simulation with a MoreauJeanOSI-Time-Stepping scheme

  see Flores/Leine/Glocker : Modeling and analysis of planar rigid multibody systems with
  translational clearance joints based on the non-smooth dynamics approach
  */

#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>

using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

#define WITH_FRICTION
// #define DISPLAY_INTER

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // parameters according to Table 1
    unsigned int nDof = 3;  // degrees of freedom for the slider crank
    double t0 = 0;          // initial computation time
    double T = 0.2;         // final computation time
    double h = 1e-5;        // time step : do not decrease, because of strong penetrations

    // geometrical characteristics
    double l1 = 0.1530;
    double l2 = 0.3060;
    double a = 0.05;
    double b = 0.025;
    double c = 0.001;

    // contact parameters
    double eN1 = 0.4;
    double eN2 = 0.4;
    double eN3 = 0.4;
    double eN4 = 0.4;
    // eN1 = 0.1;
    // eN2 = 0.1;
    // eN3 = 0.1;
    // eN4 = 0.1;
#ifdef WITH_FRICTION
    double eT1 = 0.;
    double eT2 = 0.;
    double eT3 = 0.;
    double eT4 = 0.;
    double mu1 = 0.01;
    double mu2 = 0.01;
    double mu3 = 0.01;
    double mu4 = 0.01;
#endif
    // initial conditions
    auto q0 = std::make_shared<Vector>(nDof);
    auto v0 = std::make_shared<Vector>(nDof);
    q0->zero();
    v0->zero();
    (*v0)(0) = 150.;
    (*v0)(1) = -75.;
    (*v0)(2) = -.01;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ..." << endl << endl;

    auto slider =
        std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, "SliderCrankPlugin:mass");
    slider->setComputeFGyrFunction("SliderCrankPlugin", "FGyr");
    slider->setComputeJacobianFGyrqFunction("SliderCrankPlugin", "jacobianFGyrq");
    slider->setComputeJacobianFGyrqDotFunction("SliderCrankPlugin", "jacobianFGyrqDot");
    slider->setComputeFIntFunction("SliderCrankPlugin", "FInt");
    slider->setComputeJacobianFIntqFunction("SliderCrankPlugin", "jacobianFIntq");
    slider->setComputeJacobianFIntqDotFunction("SliderCrankPlugin", "jacobianFIntqDot");

    // -------------------
    // --- Interactions---
    // -------------------
    // -- corner 1 --
#ifdef WITH_FRICTION
    auto nslaw1 =
        std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(eN1, eT1, mu1, 2);
    auto relation1 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g1", "SliderCrankPlugin:W1");
    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);

    // -- corner 2 --
    auto nslaw2 =
        std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(eN2, eT2, mu2, 2);
    auto relation2 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g2", "SliderCrankPlugin:W2");
    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation2);

    // -- corner 3 --
    auto nslaw3 =
        std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(eN3, eT3, mu3, 2);
    auto relation3 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g3", "SliderCrankPlugin:W3");
    auto inter3 = std::make_shared<siconos::modeling::Interaction>(nslaw3, relation3);

    // -- corner 4 --
    auto nslaw4 =
        std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(eN4, eT4, mu4, 2);
    auto relation4 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g4", "SliderCrankPlugin:W4");
    auto inter4 = std::make_shared<siconos::modeling::Interaction>(nslaw4, relation4);
#else
    // -- corner 1 --
    auto nslaw1 = std::make_shared<siconos::modeling::NewtonImpactNSL>(eN1);
    auto relation1 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g1", "SliderCrankPlugin:W1");
    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);

    // -- corner 2 --
    auto nslaw2 = std::make_shared<siconos::modeling::NewtonImpactNSL>(eN2);
    auto relation2 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g2", "SliderCrankPlugin:W2");
    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation2);

    // -- corner 3 --
    auto nslaw3 = std::make_shared<siconos::modeling::NewtonImpactNSL>(eN3);
    auto relation3 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g3", "SliderCrankPlugin:W3");
    auto inter3 = std::make_shared<siconos::modeling::Interaction>(nslaw3, relation3);

    // -- corner 4 --
    auto nslaw4 = std::make_shared<siconos::modeling::NewtonImpactNSL>(eN4);
    auto relation4 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g4", "SliderCrankPlugin:W4");
    auto inter4 = std::make_shared<siconos::modeling::Interaction>(nslaw4, relation4);
#endif

    // -------------
    // --- Model ---
    // -------------
    auto sliderWithClearance =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    sliderWithClearance->insertDynamicalSystem(slider);
    sliderWithClearance->link(inter1, slider);
    sliderWithClearance->link(inter2, slider);
    sliderWithClearance->link(inter3, slider);
    sliderWithClearance->link(inter4, slider);

    // ----------------
    // --- Simulation ---
    // ----------------
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(0.5);
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
#ifdef WITH_FRICTION
    auto impact = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(
        2, SICONOS_FRICTION_2D_ENUM);
#else
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>(SICONOS_LCP_LEMKE);
#endif
    impact->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    impact->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
    auto s = std::make_shared<siconos::simulation::TimeStepping>(sliderWithClearance, t);
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(200);

    auto topo = sliderWithClearance->topology();

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h) + 1;  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 27;
    Matrix dataPlot(N + 1, outputSize);

    auto q = slider->q();
    auto v = slider->velocity();

    // computation for a first consistent output
    inter1->computeOutput(t0, 0);
    inter2->computeOutput(t0, 0);
    inter3->computeOutput(t0, 0);
    inter4->computeOutput(t0, 0);

    dataPlot(0, 0) = sliderWithClearance->t0();
    dataPlot(0, 1) = (*q)(0) / (2. * M_PI);  // crank revolution
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*v)(0);
    dataPlot(0, 5) = (*v)(1);
    dataPlot(0, 6) = (*v)(2);
    dataPlot(0, 7) =
        (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) + b * cos((*q)(2)) - b) /
        c;  // y corner 1 (normalized)
    dataPlot(0, 8) =
        (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) + b * cos((*q)(2)) - b) /
        c;  // y corner 2 (normalized)
    dataPlot(0, 9) =
        (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) - b * cos((*q)(2)) + b) /
        (-c);  // y corner 3 (normalized)
    dataPlot(0, 10) =
        (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) - b * cos((*q)(2)) + b) /
        (-c);  // y corner 4 (normalized)
    dataPlot(0, 11) =
        (l1 * cos((*q)(0)) + l2 * cos((*q)(1)) - l2) / l1;          // x slider (normalized)
    dataPlot(0, 12) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1))) / c;  // y slider (normalized
    dataPlot(0, 13) = (*inter1->y(0))(0);                           // g1
    dataPlot(0, 14) = (*inter2->y(0))(0);                           // g2
    dataPlot(0, 15) = (*inter3->y(0))(0);                           // g3
    dataPlot(0, 16) = (*inter4->y(0))(0);                           // g4
    dataPlot(0, 17) = (*inter1->y(1))(0);                           // dot g1
    dataPlot(0, 18) = (*inter2->y(1))(0);                           // dot g2
    dataPlot(0, 19) = (*inter3->y(1))(0);                           // dot g3
    dataPlot(0, 20) = (*inter4->y(1))(0);                           // dot g4
    dataPlot(0, 21) = (*inter1->lambda(1))(0);                      // lambda1
    dataPlot(0, 22) = (*inter2->lambda(1))(0);                      // lambda1
    dataPlot(0, 23) = (*inter3->lambda(1))(0);                      // lambda3
    dataPlot(0, 24) = (*inter4->lambda(1))(0);                      // lambda4
    dataPlot(0, 25) = 0;
    dataPlot(0, 26) = 0;

    // --- Time loop ---
    cout << "====> Start computation ... \n";

    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();
    while ((s->hasNextEvent())) {
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"Iteration k = " << k <<std::endl;
      // std::cout <<"s->nextTime() = " <<s->nextTime()  <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;

      // std::cout << "=============== Step k ="<< k<< std::endl;
      s->advanceToEvent();
      impact->setNumericsVerboseMode(0);
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0) / (2. * M_PI);  // crank revolution
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);
      dataPlot(k, 7) =
          (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) + b * cos((*q)(2)) - b) /
          c;  // y corner 1 (normalized)
      dataPlot(k, 8) =
          (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) + b * cos((*q)(2)) - b) /
          c;  // y corner 2 (normalized)
      dataPlot(k, 9) =
          (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) - b * cos((*q)(2)) + b) /
          (c);  // y corner 3 (normalized)
      dataPlot(k, 10) =
          (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) - b * cos((*q)(2)) + b) /
          (c);  // y corner 4 (normalized)
      dataPlot(k, 11) =
          (l1 * cos((*q)(0)) + l2 * cos((*q)(1)) - l2) / l1;          // x slider (normalized)
      dataPlot(k, 12) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1))) / c;  // y slider (normalized)
      dataPlot(k, 13) = (*inter1->y(0))(0);                           // g1
      dataPlot(k, 14) = (*inter2->y(0))(0);                           // g2
      dataPlot(k, 15) = (*inter3->y(0))(0);                           // g3
      dataPlot(k, 16) = (*inter4->y(0))(0);                           // g4
      dataPlot(k, 17) = (*inter1->y(1))(0);                           // dot g1
      dataPlot(k, 18) = (*inter2->y(1))(0);                           // dot g2
      dataPlot(k, 19) = (*inter3->y(1))(0);                           // dot g3
      dataPlot(k, 20) = (*inter4->y(1))(0);                           // dot g4
      dataPlot(k, 21) = (*inter1->lambda(1))(0);                      // lambda1
      dataPlot(k, 22) = (*inter2->lambda(1))(0);                      // lambda1
      dataPlot(k, 23) = (*inter3->lambda(1))(0);                      // lambda3
      dataPlot(k, 24) = (*inter4->lambda(1))(0);                      // lambda4
      dataPlot(k, 25) = s->getNewtonNbIterations();
      auto indexSet1 = topo->indexSet(1);
      dataPlot(k, 26) = indexSet1->size();

      if (indexSet1->size() > 5) {
        impact->display();
      }
      //      if (s->nextTime() > 0.035 and (*inter1->lambda(1))(0) >0.0)
#ifdef DISPLAY_INTER
      std::cout << "=============== Step k =" << k << std::endl;
      std::cout << "Time " << s->nextTime() << std::endl;

      impact->display();
      std::cout << " (*inter1->lambda(1))(0) " << (*inter1->lambda(1))(0) << std::endl;
      std::cout << " (*inter2->lambda(1))(0) " << (*inter2->lambda(1))(0) << std::endl;
      std::cout << " (*inter3->lambda(1))(0) " << (*inter3->lambda(1))(0) << std::endl;
      std::cout << " (*inter4->lambda(1))(0) " << (*inter4->lambda(1))(0) << std::endl;
#endif

      s->processEvents();

      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("SliderCrankMoreauJeanOSI.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "SliderCrankMoreauJeanOSI.ref",
                                                      eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
