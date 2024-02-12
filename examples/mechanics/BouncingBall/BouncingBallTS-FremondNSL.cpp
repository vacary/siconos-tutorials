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

/*!\file BouncingBallTS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include "FrictionContact.hpp"
#include "TimeStepping.hpp"
#include <SiconosKernel.hpp>
#include <chrono>
using namespace std;


int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10;                  // final computation time
    double h = 0.005;                // time step
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    double R = 0.1; // Ball radius
    double m = 1; // Ball mass
    double g = 9.81; // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." <<  endl;

    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 2. / 5 * m * R * R;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0)(1) = 1.0;
    (*v0)(0) = 0.0;
    (*v0)(2) = 1.e-01;
    // -- The dynamical system --
    SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q0, v0, Mass));

    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(1) = -m * g;
    ball->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;
    double mu = 0.5;
    // Interaction ball-floor
    //
    SP::SimpleMatrix H(new SimpleMatrix(2, nDof));
    (*H)(0, 1) = 1.0;
    (*H)(1, 0) = -1.0;
    (*H)(1, 2) = -R;

    SP::NonSmoothLaw nslaw(new FremondImpactFrictionNSL(e,0.0,mu,2));
    
    SP::Relation relation(new LagrangianLinearTIR(H));

    SP::Interaction inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem bouncingBall(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->link(inter, ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));
 

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new FrictionContact(2));

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(bouncingBall, t, OSI, osnspb));

    s->setNewtonOptions(SICONOS_TS_LINEAR);
    //OSI->setGamma(3/2.0);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================


    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 12;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(1);
    SP::SiconosVector lambda = inter->lambda(1);
    SP::SiconosVector y = inter->y(0);
    SP::SiconosVector u = inter->y(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    
    dataPlot(0, 4) = (*v)(0);
    dataPlot(0, 5) = (*p)(0);
    dataPlot(0, 6) = (*y)(0);
    dataPlot(0, 7) = (*y)(1);
    dataPlot(0, 8) = (*u)(0);
    dataPlot(0, 9) = (*u)(1);
    dataPlot(0, 10) = (*lambda)(0);
    dataPlot(0, 11) = (*lambda)(1);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    while(s->hasNextEvent())
    {
      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*p)(0);
      dataPlot(k, 6) = (*y)(0);
      dataPlot(k, 7) = (*y)(1);
      dataPlot(k, 8) = (*u)(0);
      dataPlot(k, 9) = (*u)(1);
      dataPlot(k, 10) = (*lambda)(0);
      dataPlot(k, 11) = (*lambda)(1);
      
      s->nextStep();
      k++;
    }
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    cout << endl <<  "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms" << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("BouncingBallTS-FremondNSL.dat", "ascii", dataPlot, "noDim");
    ioMatrix::write("BouncingBallTS-FremondNSL.ref", "ascii", dataPlot);
    double error=0.0, eps=1e-12;
    if((error=ioMatrix::compareRefFile(dataPlot, "BouncingBallTS-FremondNSL.ref", eps)) >= 0.0
        && error > eps)
      return 1;

  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }



}
