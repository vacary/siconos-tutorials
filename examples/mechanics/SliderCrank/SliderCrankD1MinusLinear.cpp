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
  C++ input file, D1MinusLinearOSI-Time-Stepping version
  T. Schindler, V. Acary

  Slider-crank simulation with a D1MinusLinearOSI-Time-Stepping scheme

  see Flores/Leine/Glocker : Modeling and analysis of planar rigid multibody systems with
  translational clearance joints based on the non-smooth dynamics approach
  */

#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // parameters according to Table 1
    unsigned int nDof = 3;  // degrees of freedom for robot arm
    double t0 = 0.0;        // initial computation time
    double T = 0.2;         // final computation time
    // T=0.00375;

    double h = 1e-5;  // time step : do not decrease, because of strong penetrations

    // geometrical characteristics
    double l1 = 0.1530;
    double l2 = 0.3060;
    double a = 0.05;
    double b = 0.025;
    double c = 0.001;

    // contact parameters
    double e1 = 0.4;
    double e2 = 0.4;
    double e3 = 0.4;
    double e4 = 0.4;
    e1 = 0.1;
    e2 = 0.1;
    e3 = 0.1;
    e4 = 0.1;
    // double mu1 = 0.01;
    // double mu2 = 0.01;
    // double mu3 = 0.01;
    // double mu4 = 0.01;

    // initial conditions
    auto q0 = std::make_shared<Vector>(nDof);
    auto v0 = std::make_shared<Vector>(nDof);
    q0->zero();
    v0->zero();
    (*v0)(0) = 150.;
    (*v0)(1) = -75.;

    // t0 = 7e-5;
    // (*q0)(0)=  1.129178e-02;
    // (*q0)(1)= -5.777764e-03;
    // (*q0)(2)=  0.000000e+00;

    // (*v0)(0) = 1.971606e+02 ;
    // (*v0)(1) = -1.064301e+02;

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
    auto nslaw1 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e1);
    auto relation1 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g1", "SliderCrankPlugin:W1", "SliderCrankPlugin:W1dot");
    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);

    // -- corner 2 --
    auto nslaw2 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e2);
    auto relation2 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g2", "SliderCrankPlugin:W2", "SliderCrankPlugin:W2dot");
    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation2);

    // -- corner 3 --
    auto nslaw3 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e3);
    auto relation3 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g3", "SliderCrankPlugin:W3", "SliderCrankPlugin:W3dot");
    auto inter3 = std::make_shared<siconos::modeling::Interaction>(nslaw3, relation3);

    // -- corner 4 --
    auto nslaw4 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e4);
    auto relation4 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "SliderCrankPlugin:g4", "SliderCrankPlugin:W4", "SliderCrankPlugin:W4dot");
    auto inter4 = std::make_shared<siconos::modeling::Interaction>(nslaw4, relation4);

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
    auto ositype =
        siconos::integrators::D1MinusLinearOSI::Type::halfexplicit_acceleration_level;
    auto OSI = std::make_shared<siconos::integrators::D1MinusLinearOSI>(ositype);
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto force = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    auto s =
        std::make_shared<siconos::simulation::TimeSteppingD1Minus>(sliderWithClearance, t, 2);
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
    s->insertNonSmoothProblem(force, siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h) + 1;  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 35;
    Matrix dataPlot(N + 1, outputSize);

    auto q = slider->q();
    auto v = slider->velocity();

    // computation for a first consistent output
    inter1->computeOutput(t0, 0);
    inter2->computeOutput(t0, 0);
    inter3->computeOutput(t0, 0);
    inter4->computeOutput(t0, 0);

    int k = 0;
    dataPlot(k, 0) = sliderWithClearance->t0();
    dataPlot(k, 1) = (*q)(0) / (2. * M_PI);  // crank revolution
    dataPlot(k, 2) = (*q)(1);
    dataPlot(k, 3) = (*q)(2);
    dataPlot(k, 4) = (*v)(0);
    dataPlot(k, 5) = (*v)(1);
    dataPlot(k, 6) = (*v)(2);
    // std::cout << "(*q)(0)= " << (*q)(0)<< std::endl;
    // std::cout << "(*q)(1)= " << (*q)(1)<< std::endl;

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

    dataPlot(k, 13) = (*inter1->y(0))(0);       // g1
    dataPlot(k, 14) = (*inter2->y(0))(0);       // g2
    dataPlot(k, 15) = (*inter3->y(0))(0);       // g3
    dataPlot(k, 16) = (*inter4->y(0))(0);       // g4
    dataPlot(k, 17) = (*inter1->y(1))(0);       // dot g1
    dataPlot(k, 18) = (*inter2->y(1))(0);       // dot g2
    dataPlot(k, 19) = (*inter3->y(1))(0);       // dot g3
    dataPlot(k, 20) = (*inter4->y(1))(0);       // dot g4
    dataPlot(k, 21) = (*inter1->lambda(1))(0);  // lambda1
    dataPlot(k, 22) = (*inter2->lambda(1))(0);  // lambda2
    dataPlot(k, 23) = (*inter3->lambda(1))(0);  // lambda3
    dataPlot(k, 24) = (*inter4->lambda(1))(0);  // lambda4
    dataPlot(k, 25) = 0;
    dataPlot(k, 26) = 0;
    // dataPlot(k, 27) = (*inter1->lambda(2))(0) ; // lambda1_{k+1}^-
    // dataPlot(k, 28) = (*inter2->lambda(2))(0) ; // lambda1_{k+1}^-
    // dataPlot(k, 29) = (*inter3->lambda(2))(0) ; // lambda1_{k+1}^-
    // dataPlot(k, 30) = (*inter4->lambda(2))(0) ; // lambda1_{k+1}^-

    // not yet allocated
    // dataPlot(k, 31) = ( inter1->lambdaMemory(2).getSiconosVector(0) )(0); // lambda1_k^+
    // dataPlot(k, 32) = ( inter2->lambdaMemory(2).getSiconosVector(0) )(0); // lambda2_k^+
    // dataPlot(k, 33) = ( inter3->lambdaMemory(2).getSiconosVector(0) )(0); // lambda3_k^+
    // dataPlot(k, 34) = ( inter4->lambdaMemory(2).getSiconosVector(0) )(0); // lambda4_k^+

    // --- Time loop ---
    cout << "====> Start computation ... \n";

    // ==== Simulation loop - Writing without explicit event handling =====
    k++;

    auto start = std::chrono::system_clock::now();
    while ((s->hasNextEvent())) {
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"Iteration k = " << k <<std::endl;
      // std::cout <<"s->nextTime() = " <<s->nextTime()  <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;

      s->advanceToEvent();

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
      dataPlot(k, 25) = 0;
      dataPlot(k, 26) = 0;
      dataPlot(k, 27) = (*inter1->lambda(2))(0);  // lambda1_{k+1}^-
      dataPlot(k, 28) = (*inter2->lambda(2))(0);  // lambda1_{k+1}^-
      dataPlot(k, 29) = (*inter3->lambda(2))(0);  // lambda1_{k+1}^-
      dataPlot(k, 30) = (*inter4->lambda(2))(0);  // lambda1_{k+1}^-

      dataPlot(k, 31) = (inter1->lambdaMemory(2).getSiconosVector(0))(0);  // lambda1_k^+
      dataPlot(k, 32) = (inter2->lambdaMemory(2).getSiconosVector(0))(0);  // lambda2_k^+
      dataPlot(k, 33) = (inter3->lambdaMemory(2).getSiconosVector(0))(0);  // lambda3_k^+
      dataPlot(k, 34) = (inter4->lambdaMemory(2).getSiconosVector(0))(0);  // lambda4_k^+

      // std::cout << "dataPlot(k, 27)" << dataPlot(k, 27)  << std::endl;
      // std::cout << "dataPlot(k, 31)" << dataPlot(k, 31)  << std::endl;

      // std::cout <<" q->display()" <<  std::endl;
      // q->display();
      // std::cout <<" v->display()" <<  std::endl;
      // v->display();

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
    siconos::algebra::io::write("SliderCrankD1MinusLinearOSI.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "SliderCrankD1MinusLinearOSI.ref", eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
