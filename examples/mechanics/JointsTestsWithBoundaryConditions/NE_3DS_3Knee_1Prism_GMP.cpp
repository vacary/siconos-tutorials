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

/*!\file NE....cpp
  \brief \ref EMNE_MULTIBIDY - C++ input file, Time-Stepping version - O.B.

  A multibody example.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include <SolverOptions.h>

#include <KneeJointR.hpp>
#include <PrismaticJointR.hpp>
#include <SiconosKernel.hpp>
#include <chrono>

#include "GeomTools.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

#include <KneeJointR.hpp>
#include <PrismaticJointR.hpp>
#include <SiconosKernel.hpp>
#include <chrono>

#include "GeomTools.h"
using namespace std;

int main(int argc, char *argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;
    unsigned int qDim = 7;
    unsigned int nDim = 6;
    double t0 = 0;    // initial computation time
    double T = 10.0;  // final computation time
    double h = 0.01;  // time step
    int N = 1000;
    double L1 = 1.0;
    double L2 = 1.0;
    double L3 = 1.0;
    double theta = 1.0;  // theta for MoreauJeanOSI integrator
    double g = 9.81;     // Gravity
    double m = 1.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    FILE *pFile;
    pFile = fopen("data.h", "w");
    if (pFile == NULL) {
      printf("fopen exampleopen filed!\n");
      fclose(pFile);
    }

    cout << "====> Model loading ..." << endl << endl;
    // -- Initial positions and velocities --

    // First DS
    auto q10 = std::make_shared<Vector>(qDim);
    auto v10 = std::make_shared<Vector>(nDim);
    auto I1 = std::make_shared<Matrix>(3, 3);
    v10->zero();
    I1->eye();
    I1->setValue(0, 0, 0.1);
    (*q10)(0) = 0.5 * L1 / sqrt(2.0);
    (*q10)(1) = 0;
    (*q10)(2) = -0.5 * L1 / sqrt(2.0);
    double angle = M_PI / 4;
    Vector V1(3);
    V1.zero();
    V1.setValue(0, 0);
    V1.setValue(1, 1);
    V1.setValue(2, 0);
    q10->setValue(3, cos(angle / 2));
    q10->setValue(4, V1.getValue(0) * sin(angle / 2));
    q10->setValue(5, V1.getValue(1) * sin(angle / 2));
    q10->setValue(6, V1.getValue(2) * sin(angle / 2));
    // -- The dynamical system --
    auto beam1 = std::make_shared<siconos::modeling::NewtonEulerDS>(q10, v10, m, I1);
    // -- Set external forces (weight) --
    auto weight = std::make_shared<Vector>(nDof);
    (*weight)(2) = -m * g;
    beam1->setFExtPtr(weight);

    // second DS
    auto q02 = std::make_shared<Vector>(qDim);
    auto v02 = std::make_shared<Vector>(nDim);
    auto I2 = std::make_shared<Matrix>(3, 3);
    v02->zero();
    I2->eye();
    I2->setValue(0, 0, 0.1);
    (*q02)(0) = L1 / sqrt(2.0) - 0.5 * L2 / sqrt(2.0);
    (*q02)(1) = 0;
    (*q02)(2) = -L1 / sqrt(2.0) - 0.5 * L2 / sqrt(2.0);

    angle = -M_PI / 4;
    V1.zero();
    V1.setValue(0, 0);
    V1.setValue(1, 1);
    V1.setValue(2, 0);
    q02->setValue(3, cos(angle / 2));
    q02->setValue(4, V1.getValue(0) * sin(angle / 2));
    q02->setValue(5, V1.getValue(1) * sin(angle / 2));
    q02->setValue(6, V1.getValue(2) * sin(angle / 2));

    auto beam2 = std::make_shared<siconos::modeling::NewtonEulerDS>(q02, v02, m, I2);
    // -- Set external forces (weight) --
    auto weight2 = std::make_shared<Vector>(nDof);
    (*weight2)(2) = -m * g;
    beam2->setFExtPtr(weight2);

    auto q03 = std::make_shared<Vector>(qDim);
    auto v03 = std::make_shared<Vector>(nDim);
    auto I3 = std::make_shared<Matrix>(3, 3);
    v03->zero();
    I3->eye();
    I3->setValue(0, 0, 0.1);
    q03->zero();
    (*q03)(2) = -L1 * sqrt(2.0) - L1 / 2;

    angle = M_PI / 2;
    V1.zero();
    V1.setValue(0, 0);
    V1.setValue(1, 1);
    V1.setValue(2, 0);
    q03->setValue(3, cos(angle / 2));
    q03->setValue(4, V1.getValue(0) * sin(angle / 2));
    q03->setValue(5, V1.getValue(1) * sin(angle / 2));
    q03->setValue(6, V1.getValue(2) * sin(angle / 2));

    auto beam3 = std::make_shared<siconos::modeling::NewtonEulerDS>(q03, v03, m, I3);
    // -- Set external forces (weight) --
    auto weight3 = std::make_shared<Vector>(nDof);
    (*weight3)(2) = -m * g;
    beam3->setFExtPtr(weight3);
    // --------------------
    // --- Interactions ---
    // --------------------

    // Interaction with the floor
    double e = 0.9;
    auto H = std::make_shared<Matrix>(1, qDim);
    auto eR = std::make_shared<Vector>(1);
    eR->setValue(0, 2.3);
    H->zero();
    (*H)(0, 2) = 1.0;
    auto nslaw0 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation0 = std::make_shared<siconos::modeling::NewtonEulerR>();
    relation0->setJachq(H);
    relation0->setE(eR);
    cout << "main jacQH\n";
    relation0->jachq()->display();

    // --------------------
    // --- Boundary Conditions ---
    // --------------------
    auto bd = std::make_shared<siconos::modeling::BoundaryCondition>(
        siconos::modeling::BoundaryCondition::Indices{4});
    bd->setComputePrescribedVelocityFunction("Beam1Plugin", "prescribedvelocity");

    beam1->setBoundaryConditions(bd);

    // --------------------
    // --- Interactions ---
    // --------------------

    auto P = std::make_shared<Vector>(3);
    P->zero();
    auto relation1 = std::make_shared<siconos::joints::KneeJointR>(P, true, beam1);

    auto G20 = std::make_shared<Vector>(3);
    P->zero();
    P->setValue(0, L1 / 2);
    auto relation2 = std::make_shared<siconos::joints::KneeJointR>(P, false, beam1, beam2);
    P->zero();
    P->setValue(0, -L1 / 2);
    auto relation3 = std::make_shared<siconos::joints::KneeJointR>(P, false, beam2, beam3);

    auto nslaw1 = std::make_shared<siconos::modeling::EqualityConditionNSL>(
        relation1->numberOfConstraints());
    auto nslaw2 = std::make_shared<siconos::modeling::EqualityConditionNSL>(
        relation2->numberOfConstraints());
    auto nslaw3 = std::make_shared<siconos::modeling::EqualityConditionNSL>(
        relation3->numberOfConstraints());

    // Prismatic relation

    auto axe1 = std::make_shared<Vector>(3);
    axe1->zero();
    axe1->setValue(0, 1);
    auto relation4 = std::make_shared<siconos::joints::PrismaticJointR>(axe1, false, beam3);
    auto nslaw4 = std::make_shared<siconos::modeling::EqualityConditionNSL>(
        relation4->numberOfConstraints());

    // Create interations to bind relations and NSLaws.
    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);
    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation2);
    auto inter3 = std::make_shared<siconos::modeling::Interaction>(nslaw3, relation3);
    auto inter4 = std::make_shared<siconos::modeling::Interaction>(nslaw4, relation4);

    // -------------
    // --- Model ---
    // -------------
    auto myModel = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    myModel->insertDynamicalSystem(beam1);
    myModel->insertDynamicalSystem(beam2);
    myModel->insertDynamicalSystem(beam3);

    // link the interaction and the dynamical system
    myModel->link(inter1, beam1);
    myModel->link(inter2, beam1, beam2);
    myModel->link(inter3, beam2, beam3);
    myModel->link(inter4, beam3);
    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI1 = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
    auto OSI2 = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
    auto OSI3 = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::GenericMechanical>();
    osnspb->numericsSolverOptions()->dparam[0] = 1e-4;
    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(myModel, t, OSI1, osnspb);
    s->associate(OSI1, beam1);
    s->associate(OSI2, beam2);
    s->associate(OSI3, beam3);
    s->setNewtonTolerance(1e-4);
    s->setNewtonMaxIteration(50);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 15 + 7 + 6;
    Matrix dataPlot(N, outputSize);
    Matrix beam1Plot(2, 3 * N);
    Matrix beam2Plot(2, 3 * N);
    Matrix beam3Plot(2, 3 * N);

    auto q1 = beam1->q();
    auto v1 = beam1->twist();
    auto q2 = beam2->q();
    auto q3 = beam3->q();

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 0;

    fprintf(pFile, "double T[%d*%d]={", N + 1, outputSize);
    double beamTipTrajectories[6];

    auto start = std::chrono::system_clock::now();
    // N=100;
    for (k = 0; k < N; k++) {
      // solve ...
      s->advanceToEvent();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q1)(0);
      dataPlot(k, 2) = (*q1)(1);
      dataPlot(k, 3) = (*q1)(2);
      dataPlot(k, 4) = (*q1)(3);
      dataPlot(k, 5) = (*q1)(4);
      dataPlot(k, 6) = (*q1)(5);
      dataPlot(k, 7) = (*q1)(6);
      dataPlot(k, 8) = (*q2)(0);
      dataPlot(k, 9) = (*q2)(1);
      dataPlot(k, 10) = (*q2)(2);
      dataPlot(k, 11) = (*q2)(3);
      dataPlot(k, 12) = (*q2)(4);
      dataPlot(k, 13) = (*q2)(5);
      dataPlot(k, 14) = (*q2)(6);
      dataPlot(k, 15) = (*q3)(0);
      dataPlot(k, 16) = (*q3)(1);
      dataPlot(k, 17) = (*q3)(2);
      dataPlot(k, 18) = (*q3)(3);
      dataPlot(k, 19) = (*q3)(4);
      dataPlot(k, 20) = (*q3)(5);
      dataPlot(k, 21) = (*q3)(6);

      dataPlot(k, 22) = (*v1)(0);
      dataPlot(k, 23) = (*v1)(1);
      dataPlot(k, 24) = (*v1)(2);
      dataPlot(k, 25) = (*v1)(3);
      dataPlot(k, 26) = (*v1)(4);
      dataPlot(k, 27) = (*v1)(5);

      geomtools::tipTrajectories(q1, beamTipTrajectories, L1);
      beam1Plot(0, 3 * k) = beamTipTrajectories[0];
      beam1Plot(0, 3 * k + 1) = beamTipTrajectories[1];
      beam1Plot(0, 3 * k + 2) = beamTipTrajectories[2];
      beam1Plot(1, 3 * k) = beamTipTrajectories[3];
      beam1Plot(1, 3 * k + 1) = beamTipTrajectories[4];
      beam1Plot(1, 3 * k + 2) = beamTipTrajectories[5];

      geomtools::tipTrajectories(q2, beamTipTrajectories, L2);
      beam2Plot(0, 3 * k) = beamTipTrajectories[0];
      beam2Plot(0, 3 * k + 1) = beamTipTrajectories[1];
      beam2Plot(0, 3 * k + 2) = beamTipTrajectories[2];
      beam2Plot(1, 3 * k) = beamTipTrajectories[3];
      beam2Plot(1, 3 * k + 1) = beamTipTrajectories[4];
      beam2Plot(1, 3 * k + 2) = beamTipTrajectories[5];

      geomtools::tipTrajectories(q3, beamTipTrajectories, L3);
      beam3Plot(0, 3 * k) = beamTipTrajectories[0];
      beam3Plot(0, 3 * k + 1) = beamTipTrajectories[1];
      beam3Plot(0, 3 * k + 2) = beamTipTrajectories[2];
      beam3Plot(1, 3 * k) = beamTipTrajectories[3];
      beam3Plot(1, 3 * k + 1) = beamTipTrajectories[4];
      beam3Plot(1, 3 * k + 2) = beamTipTrajectories[5];

      //      printf("reaction1:%lf \n", interFloor->lambda(1)->getValue(0));

      for (unsigned int jj = 0; jj < outputSize; jj++) {
        if ((k || jj)) fprintf(pFile, ",");
        fprintf(pFile, "%f", dataPlot(k, jj));
      }
      fprintf(pFile, "\n");
      s->nextStep();
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("NE_3DS_3Knee_1Prism_GMP.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("NE_3DS_3Knee_1Prism_beam1.dat", beam1Plot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("NE_3DS_3Knee_1Prism_beam2.dat", beam2Plot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("NE_3DS_3Knee_1Prism_beam3.dat", beam3Plot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "NE_3DS_3Knee_1Prism_GMP.ref",
                                                      eps)) >= eps)
      return 1;

    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
