/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include <SiconosKernel.hpp>
#include "MechanicsFwd.hpp"
#include "BinaryCohesiveNSL.hpp"
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
    double T = 2;                  // final computation time
    double h = 1e-4;                // time step
    double position_init = 0.0;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
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
    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;

    // -- The dynamical system --
    SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q0, v0, Mass));

    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(0) = -m * g;
    /*ball->setFExtPtr(weight);*/
    ball->setComputeFExtFunction("plugins", "ballFExt");
    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;
    double sigma_c = 1e-00;
    double delta_c = 1e-04;

    // Interaction ball-floor
    //
    SP::SimpleMatrix H(new SimpleMatrix(3, nDof));
    (*H)(0, 0) = 1.0;
    (*H)(1, 1) = 1.0;
    (*H)(2, 2) = 0.0;
    
    SP::NonSmoothLaw nslaw(new BinaryCohesiveNSL(e, 0, 0, sigma_c, delta_c,3));
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
    SP::OneStepNSProblem osnspb(new CohesiveFrictionContact(3));
    osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-10; // Tolerance

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(bouncingBall, t, OSI, osnspb));

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
    SP::SiconosVector f = ball->fExt();
    SP::SiconosVector y = inter->y(0);
    SP::SiconosVector lambda = inter->lambda(1);
    int idx =0;
    dataPlot(0, idx++) = bouncingBall->t0();
    dataPlot(0, idx++) = (*q)(0);
    dataPlot(0, idx++) = (*q)(1);
    dataPlot(0, idx++) = (*v)(0);
    dataPlot(0, idx++) = (*v)(1);
    dataPlot(0, idx++) = (*p)(0);
    dataPlot(0, idx++) = (*y)(0);
    dataPlot(0, idx++) = (*y)(1);
    dataPlot(0, idx++) = (*lambda)(0);
    dataPlot(0, idx++) = (*lambda)(1);
    dataPlot(0, idx++) = 1.0;
    dataPlot(0, idx++) = (*f)(0);
    
    
    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    while(s->hasNextEvent())
    {
      std::cout << "\n\n\n next time step "<< k <<std::endl;
      s->computeOneStep();
      // --- Get values to be plotted ---
      idx=0;
      dataPlot(k, idx++) =  s->nextTime();
      dataPlot(k, idx++) = (*q)(0);
      dataPlot(k, idx++) = (*q)(1);
      dataPlot(k, idx++) = (*v)(0);
      dataPlot(k, idx++) = (*v)(1);
      dataPlot(k, idx++) = (*p)(0);
      dataPlot(k, idx++) = (*y)(0);
      dataPlot(k, idx++) = (*y)(1);
      dataPlot(k, idx++) = (*lambda)(0);
      dataPlot(k, idx++) = (*lambda)(1);
      SP::BinaryCohesiveNSL nslaw_BinaryCohesiveNSL(std::dynamic_pointer_cast<BinaryCohesiveNSL>(inter->nonSmoothLaw()));
      dataPlot(k, idx++) = nslaw_BinaryCohesiveNSL->beta(*(inter));
      dataPlot(k, idx++) = (*f)(0);
    
      std::cout << "beta = " << nslaw_BinaryCohesiveNSL->beta(*(inter)) << std::endl;
      std::cout << "f = " <<  (*f)(0) << std::endl;
      
      //getchar();
      //osnspb->display();
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
    ioMatrix::write("BouncingBallTS-BinaryCZM.dat", "ascii", dataPlot, "noDim");
    double error=0.0, eps=1e-12;
    if((error=ioMatrix::compareRefFile(dataPlot, "BouncingBallTS-BinaryCZM.ref", eps)) >= 0.0
        && error > eps)
      return 1;

  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }



}
