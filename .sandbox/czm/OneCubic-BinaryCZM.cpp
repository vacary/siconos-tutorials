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

/*!\file OneParticle-BinaryCZM.cpp
*/

#include <SiconosKernel.hpp>
#include "MechanicsFwd.hpp"
#include "BinaryCohesiveNSL.hpp"
#include <chrono>
using namespace std;

#define SYMMETRIC
//#define SHORTER_RIGHT
// #define SHORTER_LEFT


int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 0.4;                  // final computation time
    double h = 1e-4;                // time step
    //T=100*h;
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
    ball->setComputeFExtFunction("plugins", "OneCubicFExt");
    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;
    double sigma_c = 1e-00;
    double delta_c = 1e-04;
    SP::NonSmoothLaw nslaw(new BinaryCohesiveNSL(e, 0, 0, sigma_c, delta_c,3, BinaryCohesiveNSL::DOOR_SHAPE));
    //SP::NonSmoothLaw nslaw(new BinaryCohesiveNSL(e, 0, 0, sigma_c, delta_c,3, BinaryCohesiveNSL::TRIANGLE_SHAPE));
    // Interaction ball-floor
    //
    double l_x = - 1.0;
#ifdef SYMMETRIC    
    double l_y = -1. ; 
#endif
#ifdef SHORTER_RIGHT
    double l_y = -1. ;
#endif
    

    
    SP::SimpleMatrix H1(new SimpleMatrix(3, nDof));
    (*H1)(0, 0) = 1.0;
    (*H1)(0, 2) = -l_y;
    (*H1)(1, 1) = 1.0;
    (*H1)(1, 2) = l_x;
    (*H1)(2, 2) = 0.0;
    SP::Relation relation1(new LagrangianLinearTIR(H1));
    SP::Interaction inter1(new Interaction(nslaw, relation1));

    l_x = -1.0;
#ifdef SYMMETRIC 
    l_y =  1.0 ; 
#endif
#ifdef SHORTER_RIGHT
    l_y  =  .5;
#endif
    SP::SimpleMatrix H2(new SimpleMatrix(3, nDof));
    (*H2)(0, 0) = 1.0;
    (*H2)(0, 2) = -l_y;   
    (*H2)(1, 1) = 1.0;
    (*H2)(1, 2) = l_x;
    (*H2)(2, 2) = 0.0;
    SP::Relation relation2(new LagrangianLinearTIR(H2));
    SP::Interaction inter2(new Interaction(nslaw, relation2));
    
    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem bouncingBall(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->link(inter1, ball);
    bouncingBall->link(inter2, ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));
    OSI->setConstraintActivationThreshold(1e-6);

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
    unsigned int outputSize = 17;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(1);
    SP::SiconosVector f = ball->fExt();
    SP::SiconosVector y1 = inter1->y(0);
    SP::SiconosVector lambda1 = inter1->lambda(1);
    SP::SiconosVector y2 = inter2->y(0);
    SP::SiconosVector lambda2 = inter2->lambda(1);
    int idx =0;
    dataPlot(0, idx++) = bouncingBall->t0();
    dataPlot(0, idx++) = (*q)(0);
    dataPlot(0, idx++) = (*q)(1);
    dataPlot(0, idx++) = (*v)(0);
    dataPlot(0, idx++) = (*v)(1);
    dataPlot(0, idx++) = (*p)(0);
    dataPlot(0, idx++) = (*y1)(0);
    dataPlot(0, idx++) = (*y1)(1);
    dataPlot(0, idx++) = (*lambda1)(0);
    dataPlot(0, idx++) = (*lambda1)(1);
    dataPlot(0, idx++) = 1.0;
    dataPlot(0, idx++) = (*f)(0);
    dataPlot(0, idx++) = (*y2)(0);
    dataPlot(0, idx++) = (*y2)(1);
    dataPlot(0, idx++) = (*lambda2)(0);
    dataPlot(0, idx++) = (*lambda2)(1);
    dataPlot(0, idx++) = 1.0;
    
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
      dataPlot(k, idx++) = (*y1)(0);
      dataPlot(k, idx++) = (*y1)(1);
      dataPlot(k, idx++) = (*lambda1)(0);
      dataPlot(k, idx++) = (*lambda1)(1);
      SP::BinaryCohesiveNSL nslaw_BinaryCohesiveNSL(std::dynamic_pointer_cast<BinaryCohesiveNSL>(inter1->nonSmoothLaw()));
      dataPlot(k, idx++) = nslaw_BinaryCohesiveNSL->beta(*(inter1));
      dataPlot(k, idx++) = (*f)(0);
      dataPlot(k, idx++) = (*y2)(0);
      dataPlot(k, idx++) = (*y2)(1);
      dataPlot(k, idx++) = (*lambda2)(0);
      dataPlot(k, idx++) = (*lambda2)(1);
      dataPlot(k, idx++) = nslaw_BinaryCohesiveNSL->beta(*(inter2));
   
      std::cout << "beta1 = " << nslaw_BinaryCohesiveNSL->beta(*(inter1)) << std::endl;
      std::cout << "beta2 = " << nslaw_BinaryCohesiveNSL->beta(*(inter2)) << std::endl;
      std::cout << "f = " <<  (*f)(0) << std::endl;
      
      //getchar();
      //osnspb->display();
      s->nextStep();
      //getchar();
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
    ioMatrix::write("OneCubic-BinaryCZM.dat", "ascii", dataPlot, "noDim");
    double error=0.0, eps=1e-12;
#ifdef SYMMETRIC 
    if((error=ioMatrix::compareRefFile(dataPlot, "OneCubic-BinaryCZM.ref", eps)) >= 0.0
        && error > eps)
      return 1;
#endif
#ifdef SHORTER_RIGHT 
    if((error=ioMatrix::compareRefFile(dataPlot, "OneCubic-BinaryCZM-right.ref", eps)) >= 0.0
       && error > eps)
	    return 1;
#endif
  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }



}
