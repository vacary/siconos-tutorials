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

#include <SiconosKernel.hpp>
#include <chrono>

using namespace std;
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

namespace user_defined {
class MyCollisionManager : public siconos::simulation::InteractionManager {
 protected:
  /** radius of the ball */
  double _R = 0.5;

 public:
  MyCollisionManager(double R) : InteractionManager() { _R = R; }
  virtual ~MyCollisionManager() noexcept = default;

  /** Called by Simulation after updating positions prior to starting
   * the Newton loop. */
  void updateInteractions(
      std::shared_ptr<siconos::simulation::Simulation> simulation) override {
    // std::cout<< "Call to updateInteractions in MyCollisionManager" << std::endl;
    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    auto indexSet0 = simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
    for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
      auto inter = indexSet0->bundle(*ui);
      // inter->display();
      if (inter->number() == 0) {
        auto r =
            std::static_pointer_cast<siconos::modeling::Lagrangian2d2DR>(inter->relation());
        auto ds1(std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(
            indexSet0->properties(*ui).source));
        auto pc = r->pc1();
        auto nnc = r->nc();
        auto q = ds1->q();
        double angle = (*q)(2);
        // std::cout << "angle = " << angle << std::endl;
        (*pc)(0) = -_R + (*q)(0);
        (*pc)(1) = 0.0 + (*q)(1);
        (*nnc)(0) = 1.0;
        (*nnc)(1) = 0.0;
      }
    }
  }
};
}  // namespace user_defined

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;       // degrees of freedom for the ball
    double t0 = 0;               // initial computation time
    double T = 10;               // final computation time
    double h = 0.005;            // time step
    double position_init = 1.0;  // initial position for lowest bead.
    double velocity_init = 0.0;  // initial velocity for lowest bead.
    double rotation_init = 2.0;  // initial velocity for lowest bead.
    double theta = 0.5;          // theta for MoreauJeanOSI integrator
    double R = 0.5;              // Ball radius
    double m = 1;                // Ball mass
    double g = 9.81;             // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl;

    auto Mass = std::make_shared<Matrix>(nDof, nDof);
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 2. / 5 * m * R * R;

    // -- Initial positions and velocities --
    auto q0 = std::make_shared<Vector>(nDof);
    auto v0 = std::make_shared<Vector>(nDof);
    (*q0)(0) = position_init;
    (*q0)(1) = 0.0;
    (*v0)(0) = velocity_init;
    (*v0)(2) = rotation_init;

    // -- The dynamical system --
    auto ball = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, v0, Mass);

    auto q01 = std::make_shared<Vector>(nDof);
    auto v01 = std::make_shared<Vector>(nDof);
    (*q01)(0) = position_init + 2 * R + 0.1;
    (*v01)(0) = velocity_init;
    (*v01)(2) = rotation_init;

    // -- Set external forces (weight) --
    auto weight = std::make_shared<Vector>(nDof);
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(e, 0.0, 0.1, 2);

    auto relation = std::make_shared<siconos::modeling::Lagrangian2d2DR>();

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
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(2);
    // auto osnspb= std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(bouncingBall, t, OSI, osnspb);

    auto collision_manager = std::make_shared<user_defined::MyCollisionManager>(R);
    s->insertInteractionManager(collision_manager);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h) + 1;  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    Matrix dataPlot(N + 1, outputSize);

    auto q = ball->q();
    auto v = ball->velocity();
    auto p = ball->p(1);
    auto lambda = inter->lambda(1);

    dataPlot(0, 0) = s->nextTime();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*v)(2);
    dataPlot(0, 4) = (*p)(0);
    dataPlot(0, 5) = (*lambda)(0);
    dataPlot(0, 6) = (*q)(1);
    dataPlot(0, 7) = (*q)(2);
    dataPlot(0, 8) = (*v)(1);

    // --- Time loop ---
    cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    auto start = std::chrono::system_clock::now();

    while (s->hasNextEvent()) {
      // std::cout << "new time step : " << s->nextTime() <<  std::endl;

      s->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*v)(2);
      dataPlot(k, 4) = (*p)(0);
      dataPlot(k, 5) = (*lambda)(0);
      dataPlot(k, 6) = (*q)(1);
      dataPlot(k, 7) = (*q)(2);
      dataPlot(k, 8) = (*v)(1);
      s->nextStep();

      k++;
    }
    cout << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time : " << endl;
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Computation time : " << elapsed << " ms\n";

    // --- Output files ---
    cout << "====> Output file writing ...\n";
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("Ball2D_kernel_only_with_friction.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "Ball2D_kernel_only_with_friction.ref", eps)) > eps)
      return 1;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
