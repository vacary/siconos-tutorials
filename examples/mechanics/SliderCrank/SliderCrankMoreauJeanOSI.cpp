/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

  see Flores/Leine/Glocker : Modeling and analysis of planar rigid multibody
  systems with translational clearance joints based on the non-smooth dynamics
  approach
  */
#include "SCConst.h" // Parameters (geometry ...), all in namespace parameters::
#include "SiconosKernel.hpp"
#include "SolverOptions.h"
#include <chrono>

#define WITH_FRICTION
//#define DISPLAY_INTER
using namespace std;

int main(int argc, char *argv[]) {
  try {
    // ================= Creation of the model =======================

    // parameters according to Table 1
    unsigned int nDof = 3; // degrees of freedom for the slider crank
    double t0 = 0;         // initial computation time
    double T = 0.2;        // final computation time
    double h =
        1e-5; // time step : do not decrease, because of strong penetrations

    // contact parameters
    double eN1 = 0.4;
    double eN2 = 0.4;
    double eN3 = 0.4;
    double eN4 = 0.4;
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
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    q0->zero();
    v0->zero();
    (*v0)(0) = 150.;
    (*v0)(1) = -75.;
    (*v0)(2) = -.01;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ..." << endl << endl;

    SP::LagrangianDS slider(new LagrangianDS(q0, v0, "SliderCrankPlugin:mass"));
    slider->setComputeFGyrFunction("SliderCrankPlugin", "FGyr");
    slider->setComputeJacobianFGyrqFunction("SliderCrankPlugin",
                                            "jacobianFGyrq");
    slider->setComputeJacobianFGyrqDotFunction("SliderCrankPlugin",
                                               "jacobianFGyrqDot");
    slider->setComputeFIntFunction("SliderCrankPlugin", "FInt");
    slider->setComputeJacobianFIntqFunction("SliderCrankPlugin",
                                            "jacobianFIntq");
    slider->setComputeJacobianFIntqDotFunction("SliderCrankPlugin",
                                               "jacobianFIntqDot");

    // -------------------
    // --- Interactions---
    // -------------------
    // -- corner 1 --
#ifdef WITH_FRICTION
    SP::NonSmoothLaw nslaw1(new NewtonImpactFrictionNSL(eN1, eT1, mu1, 2));
    SP::Relation relation1(new LagrangianScleronomousR("SliderCrankPlugin:g1",
                                                       "SliderCrankPlugin:W1"));
    SP::Interaction inter1(new Interaction(nslaw1, relation1));

    // -- corner 2 --
    SP::NonSmoothLaw nslaw2(new NewtonImpactFrictionNSL(eN2, eT2, mu2, 2));
    SP::Relation relation2(new LagrangianScleronomousR("SliderCrankPlugin:g2",
                                                       "SliderCrankPlugin:W2"));
    SP::Interaction inter2(new Interaction(nslaw2, relation2));

    // -- corner 3 --
    SP::NonSmoothLaw nslaw3(new NewtonImpactFrictionNSL(eN3, eT3, mu3, 2));
    SP::Relation relation3(new LagrangianScleronomousR("SliderCrankPlugin:g3",
                                                       "SliderCrankPlugin:W3"));
    SP::Interaction inter3(new Interaction(nslaw3, relation3));

    // -- corner 4 --
    SP::NonSmoothLaw nslaw4(new NewtonImpactFrictionNSL(eN4, eT4, mu4, 2));
    SP::Relation relation4(new LagrangianScleronomousR("SliderCrankPlugin:g4",
                                                       "SliderCrankPlugin:W4"));
    SP::Interaction inter4(new Interaction(nslaw4, relation4));
#else
    // -- corner 1 --
    SP::NonSmoothLaw nslaw1(new NewtonImpactNSL(eN1));
    SP::Relation relation1(new LagrangianScleronomousR("SliderCrankPlugin:g1",
                                                       "SliderCrankPlugin:W1"));
    SP::Interaction inter1(new Interaction(nslaw1, relation1));

    // -- corner 2 --
    SP::NonSmoothLaw nslaw2(new NewtonImpactNSL(eN2));
    SP::Relation relation2(new LagrangianScleronomousR("SliderCrankPlugin:g2",
                                                       "SliderCrankPlugin:W2"));
    SP::Interaction inter2(new Interaction(nslaw2, relation2));

    // -- corner 3 --
    SP::NonSmoothLaw nslaw3(new NewtonImpactNSL(eN3));
    SP::Relation relation3(new LagrangianScleronomousR("SliderCrankPlugin:g3",
                                                       "SliderCrankPlugin:W3"));
    SP::Interaction inter3(new Interaction(nslaw3, relation3));

    // -- corner 4 --
    SP::NonSmoothLaw nslaw4(new NewtonImpactNSL(eN4));
    SP::Relation relation4(new LagrangianScleronomousR("SliderCrankPlugin:g4",
                                                       "SliderCrankPlugin:W4"));
    SP::Interaction inter4(new Interaction(nslaw4, relation4));
#endif

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem sliderWithClearance(
        new NonSmoothDynamicalSystem(t0, T));
    sliderWithClearance->insertDynamicalSystem(slider);
    sliderWithClearance->link(inter1, slider);
    sliderWithClearance->link(inter2, slider);
    sliderWithClearance->link(inter3, slider);
    sliderWithClearance->link(inter4, slider);

    // ----------------
    // --- Simulation ---
    // ----------------
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(0.5));
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
#ifdef WITH_FRICTION
    SP::OneStepNSProblem impact(
        new FrictionContact(2, SICONOS_FRICTION_2D_ENUM));
#else
    SP::OneStepNSProblem impact(new LCP(SICONOS_LCP_LEMKE));
#endif
    impact->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    impact->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
    SP::TimeStepping s(new TimeStepping(sliderWithClearance, t));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_TS_VELOCITY);
    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(200);

    SP::Topology topo = sliderWithClearance->topology();

    // =========================== End of model definition
    // ===========================

    // ================================= Computation
    // =================================

    int N = ceil((T - t0) / h) + 1; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 27;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = slider->q();
    SP::SiconosVector v = slider->velocity();

    // computation for a first consistent output
    inter1->computeOutput(t0, 0);
    inter2->computeOutput(t0, 0);
    inter3->computeOutput(t0, 0);
    inter4->computeOutput(t0, 0);

    int k = 0;
    dataPlot(k, 0) = sliderWithClearance->t0();
    dataPlot(k, 1) = (*q)(0) / (2. * M_PI); // crank revolution
    dataPlot(k, 2) = (*q)(1);
    dataPlot(k, 3) = (*q)(2);
    dataPlot(k, 4) = (*v)(0);
    dataPlot(k, 5) = (*v)(1);
    dataPlot(k, 6) = (*v)(2);
    // std::cout << "(*q)(0)= " << (*q)(0)<< std::endl;
    // std::cout << "(*q)(1)= " << (*q)(1)<< std::endl;

    dataPlot(k, 7) =
        (parameters::l1 * sin((*q)(0)) + parameters::l2 * sin((*q)(1)) -
         parameters::a * sin((*q)(2)) + parameters::b * cos((*q)(2)) -
         parameters::b) /
        parameters::c; // y corner 1 (normalized)
    dataPlot(k, 8) =
        (parameters::l1 * sin((*q)(0)) + parameters::l2 * sin((*q)(1)) +
         parameters::a * sin((*q)(2)) + parameters::b * cos((*q)(2)) -
         parameters::b) /
        parameters::c; // y corner 2 (normalized)
    dataPlot(k, 9) =
        (parameters::l1 * sin((*q)(0)) + parameters::l2 * sin((*q)(1)) -
         parameters::a * sin((*q)(2)) - parameters::b * cos((*q)(2)) +
         parameters::b) /
        (-parameters::c); // y corner 3 (normalized)
    dataPlot(k, 10) =
        (parameters::l1 * sin((*q)(0)) + parameters::l2 * sin((*q)(1)) +
         parameters::a * sin((*q)(2)) - parameters::b * cos((*q)(2)) +
         parameters::b) /
        (-parameters::c); // y corner 4 (normalized)

    dataPlot(k, 11) = (parameters::l1 * cos((*q)(0)) +
                       parameters::l2 * cos((*q)(1)) - parameters::l2) /
                      parameters::l1; // x slider (normalized)
    dataPlot(k, 12) =
        (parameters::l1 * sin((*q)(0)) + parameters::l2 * sin((*q)(1))) /
        parameters::c; // y slider (normalized)

    dataPlot(k, 13) = (*inter1->y(0))(0);      // g1
    dataPlot(k, 14) = (*inter2->y(0))(0);      // g2
    dataPlot(k, 15) = (*inter3->y(0))(0);      // g3
    dataPlot(k, 16) = (*inter4->y(0))(0);      // g4
    dataPlot(k, 17) = (*inter1->y(1))(0);      // dot g1
    dataPlot(k, 18) = (*inter2->y(1))(0);      // dot g2
    dataPlot(k, 19) = (*inter3->y(1))(0);      // dot g3
    dataPlot(k, 20) = (*inter4->y(1))(0);      // dot g4
    dataPlot(k, 21) = (*inter1->lambda(1))(0); // lambda1
    dataPlot(k, 22) = (*inter2->lambda(1))(0); // lambda2
    dataPlot(k, 23) = (*inter3->lambda(1))(0); // lambda3
    dataPlot(k, 24) = (*inter4->lambda(1))(0); // lambda4
    dataPlot(k, 25) = 0;
    dataPlot(k, 26) = 0;

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;

    // ==== Simulation loop - Writing without explicit event handling =====
    k = 1;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    //    while ((s->hasNextEvent()) && (k<= 3000))
    while ((s->hasNextEvent())) {

      // std::cout <<"====================================================="
      // <<std::endl; std::cout
      // <<"=====================================================" <<std::endl;
      // std::cout <<"====================================================="
      // <<std::endl; std::cout <<"Iteration k = " << k <<std::endl; std::cout
      // <<"s->nextTime() = " <<s->nextTime()  <<std::endl; std::cout
      // <<"=====================================================" <<std::endl;

      // std::cout << "=============== Step k ="<< k<< std::endl;
      s->advanceToEvent();
      impact->setNumericsVerboseMode(0);
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0) / (2. * M_PI); // crank revolution
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);
      dataPlot(k, 7) =
          (parameters::l1 * sin((*q)(0)) + parameters::l2 * sin((*q)(1)) -
           parameters::a * sin((*q)(2)) + parameters::b * cos((*q)(2)) -
           parameters::b) /
          parameters::c; // y corner 1 (normalized)
      dataPlot(k, 8) =
          (parameters::l1 * sin((*q)(0)) + parameters::l2 * sin((*q)(1)) +
           parameters::a * sin((*q)(2)) + parameters::b * cos((*q)(2)) -
           parameters::b) /
          parameters::c; // y corner 2 (normalized)
      dataPlot(k, 9) =
          (parameters::l1 * sin((*q)(0)) + parameters::l2 * sin((*q)(1)) -
           parameters::a * sin((*q)(2)) - parameters::b * cos((*q)(2)) +
           parameters::b) /
          (parameters::c); // y corner 3 (normalized)
      dataPlot(k, 10) =
          (parameters::l1 * sin((*q)(0)) + parameters::l2 * sin((*q)(1)) +
           parameters::a * sin((*q)(2)) - parameters::b * cos((*q)(2)) +
           parameters::b) /
          (parameters::c); // y corner 4 (normalized)
      dataPlot(k, 11) = (parameters::l1 * cos((*q)(0)) +
                         parameters::l2 * cos((*q)(1)) - parameters::l2) /
                        parameters::l1; // x slider (normalized)
      dataPlot(k, 12) =
          (parameters::l1 * sin((*q)(0)) + parameters::l2 * sin((*q)(1))) / parameters::c; // y slider (normalized)
      dataPlot(k, 13) = (*inter1->y(0))(0);            // g1
      dataPlot(k, 14) = (*inter2->y(0))(0);            // g2
      dataPlot(k, 15) = (*inter3->y(0))(0);            // g3
      dataPlot(k, 16) = (*inter4->y(0))(0);            // g4
      dataPlot(k, 17) = (*inter1->y(1))(0);            // dot g1
      dataPlot(k, 18) = (*inter2->y(1))(0);            // dot g2
      dataPlot(k, 19) = (*inter3->y(1))(0);            // dot g3
      dataPlot(k, 20) = (*inter4->y(1))(0);            // dot g4
      dataPlot(k, 21) = (*inter1->lambda(1))(0);       // lambda1
      dataPlot(k, 22) = (*inter2->lambda(1))(0);       // lambda1
      dataPlot(k, 23) = (*inter3->lambda(1))(0);       // lambda3
      dataPlot(k, 24) = (*inter4->lambda(1))(0);       // lambda4
      dataPlot(k, 25) = s->getNewtonNbIterations();
      SP::InteractionsGraph indexSet1 = topo->indexSet(1);
      dataPlot(k, 26) = indexSet1->size();

      if (indexSet1->size() > 5) {
        impact->display();
      }
      //      if (s->nextTime() > 0.035 and (*inter1->lambda(1))(0) >0.0)
#ifdef DISPLAY_INTER
      std::cout << "=============== Step k =" << k << std::endl;
      std::cout << "Time " << s->nextTime() << std::endl;

      impact->display();
      std::cout << " (*inter1->lambda(1))(0) " << (*inter1->lambda(1))(0)
                << std::endl;
      std::cout << " (*inter2->lambda(1))(0) " << (*inter2->lambda(1))(0)
                << std::endl;
      std::cout << " (*inter3->lambda(1))(0) " << (*inter3->lambda(1))(0)
                << std::endl;
      std::cout << " (*inter4->lambda(1))(0) " << (*inter4->lambda(1))(0)
                << std::endl;
#endif

      s->processEvents();

      k++;
    }

    cout << endl
         << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << endl;
    ;
    end = std::chrono::system_clock::now();
    int elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << "Computation time : " << elapsed << " ms" << endl;
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");

    double error = 0.0, eps = 1e-10;
    if ((error = ioMatrix::compareRefFile(
             dataPlot, "SliderCrankMoreauJeanOSI.ref", eps)) >= 0.0 &&
        error > eps)
      return 1;

  }

  catch (...) {
    Siconos::exception::process();
    return 1;
  }
}
