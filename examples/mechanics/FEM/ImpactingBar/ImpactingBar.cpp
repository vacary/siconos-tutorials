/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*!\file ImpactingBar.cpp
  V. Acary

  A Bar bouncing on the ground
  Simulation with a Time-Stepping scheme.
*/

#include <SiconosKernel.hpp>
#include <chrono>

#include "UserDefinedParameter.hpp"

using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;
namespace user =
    user_defined_ref;  // To choose the set of parameters used in the current simulation

// #define TS_PROJ
// #define TS_COMBINED

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    std::cout << "====> Model loading ...\n";

    auto ndof = user::nDof;
    double l = user::L / ndof;  // length of an element

    auto SparseMass =
        std::make_shared<Matrix>(ndof, ndof, siconos::algebra::UblasType::SPARSE, ndof);
    auto SparseStiffness =
        std::make_shared<Matrix>(ndof, ndof, siconos::algebra::UblasType::SPARSE, 3 * ndof);

    SparseMass->setValue(0, 0, 1.0 / 3.0);
    SparseMass->setValue(0, 1, 1.0 / 6.0);
    SparseStiffness->setValue(0, 0, 1.0);
    SparseStiffness->setValue(0, 1, -1.0);

    for (unsigned int i = 1; i < ndof - 1; i++) {
      SparseMass->setValue(i, i, 2.0 / 3.0);
      SparseMass->setValue(i, i - 1, 1.0 / 6.0);
      SparseMass->setValue(i, i + 1, 1.0 / 6.0);

      SparseStiffness->setValue(i, i, 2.0);
      SparseStiffness->setValue(i, i - 1, -1.0);
      SparseStiffness->setValue(i, i + 1, -1.0);
    }

    SparseMass->setValue(ndof - 1, ndof - 1, 1.0 / 3.0);
    SparseMass->setValue(ndof - 1, ndof - 2, 1.0 / 6.0);

    SparseStiffness->setValue(ndof - 1, ndof - 2, -1.0);
    SparseStiffness->setValue(ndof - 1, ndof - 1, 1.0);

    *SparseMass *= user::rho * user::S * l;
    *SparseStiffness *= user::E * user::S / l;

    //      SparseMass->display();
    //      SparseStiffness->display();

    // -- Initial positions and velocities --
    auto q0 = std::make_shared<Vector>(ndof, user::position_init);
    auto v0 = std::make_shared<Vector>(ndof, user::velocity_init);

    // -- The dynamical system --
    auto bar = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, v0, SparseMass);

    // -- Set stiffness matrix (weight) --
    bar->setKPtr(SparseStiffness);

    // -- Set external forces (weight) --
    // auto weight= std::make_shared<Vector>(ndof,-g*rho*S/l);
    auto weight = std::make_shared<Vector>(ndof, 0.0);
    bar->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.0;

    // Interaction bar-floor
    //
    auto H = std::make_shared<Matrix>(1, ndof);
    (*H)(0, 0) = 1.0;

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H);

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // -------------
    // --- Model ---
    // -------------
    auto impactingBar =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(user::t0, user::T);

    // add the dynamical system in the non smooth dynamical system
    impactingBar->insertDynamicalSystem(bar);

    // link the interaction and the dynamical system
    impactingBar->link(inter, bar);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
#ifdef TS_PROJ
    auto OSI =
        std::make_shared<siconos::integrators::MoreauJeanDirectProjectionOSI>(user::theta);
    OSI->setDeactivateYPosThreshold(1e-05);
    OSI->setDeactivateYVelThreshold(0.0);
    OSI->setActivateYPosThreshold(1e-09);
    OSI->setActivateYVelThreshold(100.0);
#else
#ifdef TS_COMBINED
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanCombinedProjectionOSI>(theta);
#else
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(user::theta, 0.5);
#endif
#endif
    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(user::t0, user::h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
#ifdef TS_PROJ
    auto position =
        std::make_shared<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>();
    auto s = std::make_shared<siconos::simulation::TimeStepping> DirectProjection(
        impactingBar, t, OSI, osnspb, position, 0);
    s->setProjectionMaxIteration(10);
    s->setConstraintTolUnilateral(1e-10);
    s->setConstraintTol(1e-10);

#else
#ifdef TS_COMBINED
    auto position =
        std::make_shared<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>(
            SICONOS_MLCP_ENUM);
    SP::TimeSteppingCombinedProjection s =
        std::make_shared<siconos::simulation::TimeStepping> CombinedProjection(
            impactingBar, t, OSI, osnspb, position, 2);
    s->setProjectionMaxIteration(500);
    s->setConstraintTolUnilateral(1e-10);
    s->setConstraintTol(1e-10);
#else
    auto s = std::make_shared<siconos::simulation::TimeStepping>(impactingBar, t, OSI, osnspb);
#endif
#endif

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    int N = floor((user::T - user::t0) / user::h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 12;
    Matrix dataPlot(N, outputSize);

    auto q = bar->q();
    auto v = bar->velocity();
    auto p = bar->p(1);
    auto lambda = inter->lambda(1);

    dataPlot(0, 0) = impactingBar->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 7) = (*q)(ndof - 1);
    dataPlot(0, 8) = (*v)(ndof - 1);
    dataPlot(0, 9) = (*q)((ndof) / 2);
    dataPlot(0, 10) = (*v)((ndof) / 2);

    auto tmp = std::make_shared<Vector>(ndof);

    prod(*SparseStiffness, *q, *tmp, true);
    double potentialEnergy = inner_prod(*q, *tmp);
    prod(*SparseMass, *v, *tmp, true);
    double kineticEnergy = inner_prod(*v, *tmp);
    double impactEnergy = 0.0;

    dataPlot(0, 5) = potentialEnergy;
    dataPlot(0, 6) = kineticEnergy;
    dataPlot(0, 11) = impactEnergy;

    //    std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
    //     std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();
    //    while (s->nextTime() < T)
    while (k < N) {
      s->computeOneStep();

      //       std::cout << "position"  << std::endl;
      //       q->display();
      //       std::cout << "velocity"  << std::endl;
      //       v->display();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0) / user::h;
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 7) = (*q)(ndof - 1);
      dataPlot(k, 8) = (*v)(ndof - 1);
      dataPlot(k, 9) = (*q)((ndof) / 2);
      dataPlot(k, 10) = (*v)((ndof) / 2);

      prod(*SparseStiffness, *q, *tmp, true);
      potentialEnergy = inner_prod(*q, *tmp);
      prod(*SparseMass, *v, *tmp, true);
      kineticEnergy = inner_prod(*v, *tmp);

      dataPlot(k, 5) = potentialEnergy;
      dataPlot(k, 6) = kineticEnergy;

      double v_p_theta = user::theta * (*v)(0) + (1 - user::theta) * dataPlot(k - 1, 2);
      impactEnergy = v_p_theta * (*lambda)(0);
      dataPlot(k, 11) = impactEnergy;

      //      std::cout << "q" << std::endl;
      //       q->display();

      //       std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
      //       std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;

      s->nextStep();

      k++;
    }
    double totalImpactEnergy = 0.0;
    for (int p = 0; p < k; p++) {
      totalImpactEnergy = totalImpactEnergy + dataPlot(p, 11);
    }
    printf("totalImpactEnergy = %e", totalImpactEnergy);

    auto end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("ImpactingBar.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-11;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "ImpactingBar.ref", eps)) >=
        eps)
      return 1;

    return 0;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
