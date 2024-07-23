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

/*!\file BouncingBallNETS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, O. Bonnefon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include <SiconosKernel.hpp>
#include <chrono>

using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

#define WITH_PROJ
#define WITH_FC3D
using namespace std;
#ifdef WITH_FC3D
#define R_CLASS NewtonEuler3DR
#else
#define R_CLASS NewtonEuler1DR
#endif

namespace user_defined {
class my_NewtonEulerR : public siconos::modeling::R_CLASS {
  double _sBallRadius;

 public:
  my_NewtonEulerR(double radius) : R_CLASS(), _sBallRadius(radius) {};

  virtual void computeOutput(double t, siconos::modeling::Interaction& inter,
                             unsigned int derivativeNumber) override {
    auto& DSlink = inter.linkToDSVariables();
    if (derivativeNumber == 0) {
      computeh(t, *DSlink[NewtonEulerR::q0], *inter.y(0));
    } else {
      R_CLASS::computeOutput(t, inter, derivativeNumber);
    }
  }

  void computeh(double time, const siconos::algebra::BlockVector& q0,
                siconos::algebra::SiconosVector& y) override {
    std::cout << "my_NewtonEulerR:: computeh \n";
    std::cout << "q0.size() = " << q0.size() << "\n";
    double height = q0.getValue(0) - _sBallRadius - q0.getValue(7);
    // std::cout <<"my_NewtonEulerR:: computeh _jachq" << std:: endl;
    // _jachq->display();
    y.setValue(0, height);
    _Nc->setValue(0, 1);
    _Nc->setValue(1, 0);
    _Nc->setValue(2, 0);
    _Pc1->setValue(0, q0.getValue(0) - _sBallRadius);
    _Pc1->setValue(1, q0.getValue(1));
    _Pc1->setValue(2, q0.getValue(2));

    _Pc2->setValue(0, q0.getValue(7));
    _Pc2->setValue(1, q0.getValue(8));
    _Pc2->setValue(2, q0.getValue(9));
    // printf("my_NewtonEulerR N, Pc\n");
    _Nc->display();
    _Pc1->display();
    _Pc2->display();
    std::cout << "my_NewtonEulerR:: computeh ends" << std::endl;
  }
};
}  // namespace user_defined

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;       // degrees of freedom for the ball
    unsigned int qDim = 7;       // degrees of freedom for the ball
    unsigned int nDim = 6;       // degrees of freedom for the ball
    double t0 = 0;               // initial computation time
    double T = 10.0;             // final computation time
    double h = 0.001;            // time step
    double position_init = 1.0;  // initial position
    double velocity_init = 0.0;  // initial velocity
    double omega_initx = 0.0;    // initial angular velocity
    double omega_initz = 0.0;    // initial angular velocity
    double theta = 1.0;          // theta for MoreauJeanOSI integrator
    double m = 1;                // Ball mass
    double g = 10.0;             // Gravity
    double radius = 0.1;
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    std::cout << "====> Model loading ...\n";

    // -- Initial positions and velocities --
    auto q0 = std::make_shared<Vector>(qDim);
    auto v0 = std::make_shared<Vector>(nDim);
    auto I = std::make_shared<Matrix>(3, 3);
    v0->zero();
    q0->zero();

    I->eye();
    (*q0)(0) = position_init;
    /*initial quaternion equal to (1,0,0,0)*/
    (*q0)(3) = 1.0;

    (*v0)(0) = velocity_init;
    (*v0)(3) = omega_initx;
    (*v0)(5) = omega_initz;
    // -- The dynamical system --
    auto ball = std::make_shared<siconos::modeling::NewtonEulerDS>(q0, v0, m, I);

    // -- Set external forces (weight) --
    auto weight = std::make_shared<Vector>(nDof);
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

    // siconos::modeling::BoundaryCondition::Indices bdindex = {0, 3, 5};
    auto bd = std::make_shared<siconos::modeling::BoundaryCondition>(
        siconos::modeling::BoundaryCondition::Indices{0, 3, 5});
    bd->setComputePrescribedVelocityFunction("BallOnMovingPlanePlugin", "prescribedvelocity3");
    ball->setBoundaryConditions(bd);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

#ifdef WITH_FC3D
    auto nslaw0 = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(e, e, 0.6, 3);
#else
    auto nslaw0 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
#endif

    auto relation0 = std::make_shared<user_defined::my_NewtonEulerR>(radius);
    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw0, relation0);

    // -------------
    // --- Model ---
    // -------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
#ifdef WITH_PROJ
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanDirectProjectionOSI>(theta);
#else
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
#endif
    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::GenericMechanical>();
#ifdef WITH_PROJ
    auto osnspb_pos =
        std::make_shared<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>(
            SICONOS_MLCP_ENUM, 1.0);
#endif
    // -- (4) Simulation setup with (1) (2) (3)
#ifdef WITH_PROJ
    auto s = std::make_shared<siconos::simulation::TimeSteppingDirectProjection>(
        bouncingBall, t, OSI, osnspb, osnspb_pos);
    s->setProjectionMaxIteration(20);
    s->setConstraintTolUnilateral(1e-08);
    s->setConstraintTol(1e-08);
#else
    auto s = std::make_shared<siconos::simulation::TimeStepping>(bouncingBall, t, OSI, osnspb);
#endif
    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(10);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 20;
    Matrix dataPlot(N + 1, outputSize);

    auto q = ball->q();
    auto v = ball->velocity();
    auto p = ball->p(1);

    // auto lambda = inter->lambda(1);
    // auto y = inter->y(0);

    auto reaction = ball->reactionToBoundaryConditions();

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*reaction)(0);
    dataPlot(0, 5) = acos((*q)(3));
    // dataPlot(0, 6) = relation0->contactForce()->norm2();

    dataPlot(0, 7) = (*q)(0);
    dataPlot(0, 8) = (*q)(1);
    dataPlot(0, 9) = (*q)(2);
    dataPlot(0, 10) = (*q)(3);
    dataPlot(0, 11) = (*q)(4);
    dataPlot(0, 12) = (*q)(5);
    dataPlot(0, 13) = (*q)(6);

    dataPlot(0, 14) = (*v)(0);
    dataPlot(0, 15) = (*v)(1);
    dataPlot(0, 16) = (*v)(2);
    dataPlot(0, 17) = (*v)(3);
    dataPlot(0, 18) = (*v)(4);
    dataPlot(0, 19) = (*v)(5);

    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      // std::cout << "step " << k << std::endl;
      //        s->computeOneStep();
      s->advanceToEvent();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*reaction)(0);

      /* fix issue with machine accuracy */
      if ((*q)(3) > 1.0) {
        dataPlot(k, 5) = acos(1.0);
      } else if ((*q)(3) < -1.0) {
        dataPlot(k, 5) = acos(-1.0);
      } else {
        dataPlot(k, 5) = acos((*q)(3));
      }
      // dataPlot(k, 6) = relation0->contactForce()->norm2();
      dataPlot(k, 7) = (*q)(0);
      dataPlot(k, 8) = (*q)(1);
      dataPlot(k, 9) = (*q)(2);
      dataPlot(k, 10) = (*q)(3);
      dataPlot(k, 11) = (*q)(4);
      dataPlot(k, 12) = (*q)(5);
      dataPlot(k, 13) = (*q)(6);

      dataPlot(k, 14) = (*v)(0);
      dataPlot(k, 15) = (*v)(1);
      dataPlot(k, 16) = (*v)(2);
      dataPlot(k, 17) = (*v)(3);
      dataPlot(k, 18) = (*v)(4);
      dataPlot(k, 19) = (*v)(5);

      s->nextStep();

      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("BallNewtonEuler.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // Comparison with a reference file
    cout << "====> Comparison with a reference file ...\n";
#ifdef WITH_PROJ
    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "BallNewtonEuler-WITHPROJ.ref",
                                                      eps)) > eps)
      return 1;
#else
    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "BallNewtonEuler.ref", eps)) >
        eps)
      return 1;
#endif

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
