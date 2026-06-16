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

#include <BinaryCohesiveNSL.hpp>
#include <SiconosKernel.hpp>
//#include "BinaryCohesiveNSL.hpp"
#include <chrono>
using namespace std;

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

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

    Matrix mass{nDof, nDof};
    mass.setZero();
    mass(0, 0) = m;
    mass(1, 1) = m;
    mass(2, 2) = 2. / 5 * m * R * R;

    // -- Initial positions and velocities --
    Vector q0{nDof};
    q0.setZero();
    q0(0) = position_init;
    Vector v0{nDof};
    v0.setZero();
    v0(0) = velocity_init;

    // -- The dynamical system --
    auto ball = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
        q0, v0, mass, siconos::algebra::alias_t);
    // -- Set external forces (weight) --
    Vector weight{nDof};
    weight.setZero();
    weight(0) = -m * g;
    ball->setConstantFext(weight, siconos::algebra::alias_t);
    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;
    double sigma_c = 1e-00;
    double delta_c = 1e-04;

    // Interaction ball-floor
    //
    Matrix H{1, nDof};
    H.setZero();
    H(0, 0) = 1.0;
    

    auto nslaw = std::make_shared<siconos::mechanics::czm::BinaryCohesiveNSL>(e, 0, 0, sigma_c, delta_c,3);
    auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H);
    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // -------------
    // --- Model ---
    // -------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->link(inter, ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);


    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();    
    osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-10; // Tolerance

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(bouncingBall, t, OSI, osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================


    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 12;
    Matrix dataPlot(N + 1, outputSize);

    auto q = ball->q_read();
    auto v = ball->velocity_read();
    auto p = *(ball->p(1));
    auto f = ball->totalForces(); //fext()?
    auto y = *(inter->y(0));
    auto lambda = *(inter->lambda(1));
    
    int idx =0;
    dataPlot(0, idx++) = bouncingBall->t0();
    dataPlot(0, idx++) = q(0);
    dataPlot(0, idx++) = q(1);
    dataPlot(0, idx++) = v(0);
    dataPlot(0, idx++) = v(1);
    dataPlot(0, idx++) = p(0);
    dataPlot(0, idx++) = y(0);
    dataPlot(0, idx++) = y(1);
    dataPlot(0, idx++) = lambda(0);
    dataPlot(0, idx++) = lambda(1);
    dataPlot(0, idx++) = 1.0;
    dataPlot(0, idx++) = f(0);
    
    
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
      dataPlot(k, idx++) = q(0);
      dataPlot(k, idx++) = q(1);
      dataPlot(k, idx++) = v(0);
      dataPlot(k, idx++) = v(1);
      dataPlot(k, idx++) = p(0);
      dataPlot(k, idx++) = y(0);
      dataPlot(k, idx++) = y(1);
      dataPlot(k, idx++) = lambda(0);
      dataPlot(k, idx++) = lambda(1);
      auto nslaw_BinaryCohesiveNSL(std::dynamic_pointer_cast<siconos::mechanics::czm::BinaryCohesiveNSL>(inter->nonSmoothLaw()));
      dataPlot(k, idx++) = nslaw_BinaryCohesiveNSL->beta(*(inter));
      dataPlot(k, idx++) = f(0);
    
      std::cout << "beta = " << nslaw_BinaryCohesiveNSL->beta(*(inter)) << std::endl;
      std::cout << "f = " <<  f(0) << std::endl;
      
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
    siconos::algebra::io::write("BouncingBallTS-BinaryCZM.dat",  dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error=0.0, eps=1e-12;
    if((error= siconos::algebra::io::compareRefFile(dataPlot, "BouncingBallTS-BinaryCZM.ref", eps)) >= 0.0
        && error > eps)
      return 1;

  }

  catch(...)
  {
    siconos::exception::process();
    return 1;
  }



}
