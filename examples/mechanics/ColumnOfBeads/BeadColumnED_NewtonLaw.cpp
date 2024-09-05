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

/*!\file BouncingBallED.cpp
  \brief \ref EMBouncingBall - C++ input file, Event-Driven version - V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with an Event-Driven scheme.
*/

#include <SiconosKernel.hpp>
#include <chrono>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

constexpr double PI = 3.14159265;
constexpr double g = 9.81;  // Gravity

using namespace std;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDofBall = 1;  // degrees of freedom of ball 1
    double Height = 0.05;       // Distance between impactor balls and monodisperse balls
    // Balls in the tapered chain
    unsigned int NumberBalls = 10;  // Number
    double R_base_taper = 0.005;    // Base radii of the tapered chain
    double q_taper = 0.0;           // Tapering factor
    // Material properties of balls
    double mass_density = 7780;  // mass density
    double Res_BallBall = 1.0;   // Restitution coefficient at contacts ball-ball
    double Res_BallWall = 1.0;   // Restitution coefficient at contacts ball-wall
    double PowCompLaw = 1.5;     // Power of the compliance law: 1.0 for linear contact and 3/2
                                 // for the Hertzian contact
    string TypeContactLaw = "BiStiffness";  // Type of compliance contact law
    // Parameters for the global simulation
    double t0 = 0;                  // initial computation time
    double T = 0.5;                 // final computation time
    double h = 0.001;               // time step
    unsigned int Npointsave = 500;  // Number of data points to be saved
    //---------------------------------------
    // ---- Configuration of chaines
    //--------------------------------------
    //************* Balls ******************
    double NumberContacts = NumberBalls;  // Number of contacts
    //(1) Radius of balls
    auto RadiusBalls = std::make_shared<Vector>(NumberBalls);
    for (unsigned int k = 0; k < NumberBalls; ++k) {
      (*RadiusBalls)(k) = (pow(double(1.0 - q_taper), int(k + 1))) * R_base_taper;
    }
    // (2) Mass of balls
    auto MassBalls = std::make_shared<Vector>(NumberBalls);
    for (unsigned int id = 0; id < NumberBalls; ++id) {
      (*MassBalls)(id) = (4.0 / 3.0) * PI * pow((*RadiusBalls)(id), 3) * mass_density;
    }
    // (3) Initial position of balls
    // For the impactor balls
    auto InitPosBalls = std::make_shared<Vector>(NumberBalls);
    (*InitPosBalls)(0) = Height + (*RadiusBalls)(0);
    for (unsigned int j = 1; j < NumberBalls; ++j) {
      (*InitPosBalls)(j) = (*InitPosBalls)(j - 1) + (*RadiusBalls)(j - 1) + (*RadiusBalls)(j);
    }
    // (4) Initial velocity of balls
    auto InitVelBalls = std::make_shared<Vector>(NumberBalls);
    for (unsigned int i = 0; i < NumberBalls; ++i) {
      (*InitVelBalls)(i) = 0.0;
    }
    //****************** Contacts ******************
    // (1) Restitution coefficient at contacts
    auto ResCofContacts = std::make_shared<Vector>(NumberContacts);
    auto ElasCofContacts = std::make_shared<Vector>(NumberContacts);
    for (unsigned int id = 0; id < NumberContacts; ++id) {
      if (id == 0)  // contact ball-wall
      {
        (*ResCofContacts)(id) = Res_BallWall;
      } else  // contact ball-ball
      {
        (*ResCofContacts)(id) = Res_BallBall;
      }
    }

    // // Display and save the configuration of the chain simulated
    // cout << "Configuation of ball chains\n";
    // cout.precision(15);
    // cout << "Radius of balls: \n";
    // RadiusBalls->display();
    // cout << "Mass of balls: \n";
    // MassBalls->display();
    // cout << "Initial position of balls: \n";
    // InitPosBalls->display();
    // cout << "Initial velocity of balls: \n";
    // InitVelBalls->display();
    // cout << "Restitution coefficient at contacts:\n";
    // ResCofContacts->display();
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ..." << endl << endl;
    std::vector<std::shared_ptr<siconos::modeling::LagrangianLinearTIDS>> VecOfallDS;
    double _Rball, _massBall, _Pos0Ball, _Vel0Ball;
    // -------------
    // --- Model ---
    // -------------
    auto BallChain = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // ----------------
    // --- Simulation ---
    // ----------------
    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::LsodarOSI>();

    for (unsigned int i = 0; i < NumberBalls; ++i) {
      _Rball = (*RadiusBalls)(i);      // radius of the ball
      _massBall = (*MassBalls)(i);     // mass of the ball
      _Pos0Ball = (*InitPosBalls)(i);  // initial position of the ball
      _Vel0Ball = (*InitVelBalls)(i);  // initial velocity of the ball
      // Declaration of the DS in Siconos
      auto MassBall = std::make_shared<Matrix>(nDofBall, nDofBall);
      (*MassBall)(0, 0) = _massBall;
      // -- Initial positions and velocities --
      auto q0Ball = std::make_shared<Vector>(nDofBall);
      auto v0Ball = std::make_shared<Vector>(nDofBall);
      (*q0Ball)(0) = _Pos0Ball;
      (*v0Ball)(0) = _Vel0Ball;
      // -- The dynamical system --
      auto ball =
          std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0Ball, v0Ball, MassBall);
      // -- Set external forces (weight1) --
      auto FextBall = std::make_shared<Vector>(nDofBall);
      (*FextBall)(0) = -_massBall * g;
      ball->setFExtPtr(FextBall);
      //
      VecOfallDS.push_back(ball);
      BallChain->insertDynamicalSystem(ball);
    }
    // --------------------
    // --- Interactions ---
    // --------------------
    double ResCoef, Stiff, ElasPow;
    for (unsigned int j = 0; j < NumberContacts; ++j) {
      ResCoef = (*ResCofContacts)(j);

      std::shared_ptr<Vector> E;
      std::shared_ptr<Matrix> H;
      if (j == 0)  // for contact wall-ball
      {
        H = std::make_shared<Matrix>(1, nDofBall);
        (*H)(0, 0) = 1.0;
        E = std::make_shared<Vector>(1);
        (*E)(0) = -(*RadiusBalls)(j);

      } else  // For ball-ball contact
      {
        H = std::make_shared<Matrix>(1, (nDofBall + nDofBall));
        (*H)(0, 0) = -1.0;
        (*H)(0, 1) = 1.0;
        E = std::make_shared<Vector>(1);
        (*E)(0) = -1.0 * ((*RadiusBalls)(j - 1) + (*RadiusBalls)(j));
      }
      //
      auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(ResCoef);
      auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H, E);
      auto interaction = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);
      if (j == 0)  // for contact wall-ball
        BallChain->link(interaction, VecOfallDS[j]);
      else  // For ball-ball contact
        BallChain->link(interaction, VecOfallDS[j - 1], VecOfallDS[j]);
    }

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // -- (3) Non smooth problem --
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto acceleration = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::EventDriven>(BallChain, t);
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_ED_IMPACT);
    s->insertNonSmoothProblem(acceleration, siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC);

    // =========================== End of model definition ===========================
    //----------------------------------- Initialization-------------------------------

    auto DSG0 = BallChain->topology()->dSG(0);
    auto IndexSet0 = BallChain->topology()->indexSet(0);

    // // Display topology of the system
    // cout << "Number of vectices of IndexSet0: " << IndexSet0->size() << endl;
    // cout << "Number of vectices of DSG0: " << DSG0->size() << endl;
    //
    auto eventsManager = s->eventsManager();
    // ================================= Computation =================================
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 2 * NumberBalls + 1;
    Matrix dataPlot(Npointsave, outputSize);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    bool nonSmooth = false;
    unsigned int NumberOfEvents = 0;
    unsigned int NumberOfNSEvents = 0;
    unsigned int k = 0;
    siconos::graphs::DynamicalSystemsGraph::VIterator ui, uiend;
    auto start = std::chrono::system_clock::now();
    //====================================================================
    while ((k < Npointsave) & (s->hasNextEvent())) {
      dataPlot(k, 0) = s->startingTime();
      // Save state of the balls
      unsigned int col_pos = 1;
      unsigned int col_vel = NumberBalls + 1;
      for (boost::tie(ui, uiend) = DSG0->vertices(); ui != uiend; ++ui) {
        auto ds = DSG0->bundle(*ui);
        auto lag_ds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds);
        auto q = lag_ds->q();
        auto v = lag_ds->velocity();
        dataPlot(k, col_pos) = (*q)(0);
        dataPlot(k, col_vel) = (*v)(0);
        col_pos++;
        col_vel++;
      }
      ++k;
      s->advanceToEvent();  // run simulation from one event to the next
      if (eventsManager->nextEvent()->getType() == siconos::simulation::EventType::NS) {
        nonSmooth = true;
      };
      //
      s->processEvents();  // process events
      if (nonSmooth) {
        dataPlot(k, 0) = s->startingTime();
        // Save state of the balls
        unsigned int col_pos = 1;
        unsigned int col_vel = NumberBalls + 1;
        for (boost::tie(ui, uiend) = DSG0->vertices(); ui != uiend; ++ui) {
          auto ds = DSG0->bundle(*ui);
          auto lag_ds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds);
          const auto& q = lag_ds->qMemory().getSiconosVector(1);
          const auto& v = lag_ds->velocityMemory().getSiconosVector(1);
          dataPlot(k, col_pos) = q(0);
          dataPlot(k, col_vel) = v(0);
          col_pos++;
          col_vel++;
        }
        nonSmooth = false;
        ++NumberOfNSEvents;
        ++NumberOfEvents;

        ++k;
      }
      // --- Get values to be plotted ---
      ++NumberOfEvents;
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("BeadColumnED_NewtonLaw.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "BeadColumnED_NewtonLaw.ref",
                                                      eps)) > eps)
      return 1;
    return 0;
  } catch (...) {
    siconos::exception::process();
    return 1;
  }
}
