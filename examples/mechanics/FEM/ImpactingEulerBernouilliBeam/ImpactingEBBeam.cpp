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

/*!\file ImpactingBar.cpp
  V. Acary

  A Bar bouncing on the ground
  Simulation with a Time-Stepping scheme.
*/

#include <chrono>

#include "SiconosAlgebraProd.hpp"
#include "SiconosKernel.hpp"

using namespace std;

static void addElementaryStiffnessMatrix(
    std::shared_ptr<siconos::algebra::SiconosMatrix> SparseStiffness, int elementNumber,
    int nDof, double elementLength) {
  int dofStart = elementNumber * 2;
  int ndof_per_element = 4;

  std::shared_ptr<siconos::algebra::SiconosMatrix> Ke(new SiconosMatrix(4, 4));
  (*Ke)(0, 0) = 12.;
  (*Ke)(0, 1) = 6. * elementLength;
  (*Ke)(0, 2) = -12.;
  (*Ke)(0, 3) = 6. * elementLength;

  (*Ke)(1, 1) = 4. * elementLength * elementLength;
  (*Ke)(1, 2) = -6. * elementLength;
  (*Ke)(1, 3) = 2. * elementLength * elementLength;

  (*Ke)(2, 2) = 12.;
  (*Ke)(2, 3) = -6. * elementLength;

  (*Ke)(3, 3) = 4. * elementLength * elementLength;

  for (unsigned int i = 0; i < 4; i++) {
    for (unsigned int j = 1; j < 4; j++) (*Ke)(j, i) = (*Ke)(i, j);
  }

  // Ke->display();
  std::shared_ptr<siconos::algebra::SiconosMatrix> ElementStiffness(
      new SiconosMatrix(nDof, nDof, siconos::SPARSE, 10 * nDof));
  for (unsigned int i = 0; i < ndof_per_element; i++) {
    for (unsigned int j = 0; j < ndof_per_element; j++)
      ElementStiffness->setValue(i + dofStart, j + dofStart, (*Ke)(i, j));
  }
  // ElementStiffness->display();

  (*SparseStiffness) += (*ElementStiffness);
}

static void addElementaryMassMatrix(
    std::shared_ptr<siconos::algebra::SiconosMatrix> SparseMass, int elementNumber, int nDof,
    double elementLength, bool lumpedMass) {
  int dofStart = elementNumber * 2;
  int ndof_per_element = 4;

  std::shared_ptr<siconos::algebra::SiconosMatrix> Me(new SiconosMatrix(4, 4));

  if (!lumpedMass) {
    (*Me)(0, 0) = 156.;
    (*Me)(0, 1) = 22. * elementLength;
    (*Me)(0, 2) = 54.;
    (*Me)(0, 3) = -13. * elementLength;

    (*Me)(1, 1) = 4. * elementLength * elementLength;
    (*Me)(1, 2) = 13. * elementLength;
    (*Me)(1, 3) = -3. * elementLength * elementLength;

    (*Me)(2, 2) = 156.;
    (*Me)(2, 3) = -22. * elementLength;

    (*Me)(3, 3) = 4. * elementLength * elementLength;

    for (unsigned int i = 0; i < 4; i++) {
      for (unsigned int j = 1; j < 4; j++) (*Me)(j, i) = (*Me)(i, j);
    }
  } else {
    double alpha = 1 / 50.;
    (*Me)(0, 0) = 0.5 * 420.;
    (*Me)(1, 1) = alpha * elementLength * elementLength * 420.;
    (*Me)(2, 2) = 0.5 * 420.;
    (*Me)(3, 3) = alpha * elementLength * elementLength * 420.;
  }

  // Me->display();

  std::shared_ptr<siconos::algebra::SiconosMatrix> ElementMass(
      new SiconosMatrix(nDof, nDof, siconos::SPARSE, 6 * nDof));
  for (unsigned int i = 0; i < ndof_per_element; i++) {
    for (unsigned int j = 0; j < ndof_per_element; j++)
      ElementMass->setValue(i + dofStart, j + dofStart, (*Me)(i, j));
  }
  // ElementStiffness->display();

  (*SparseMass) += (*ElementMass);
}

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
#include "UserDefinedParameter.hpp"

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;

    double l = L / nElement;  // length of an element

    cout << "beam length = " << l << endl;

    unsigned int nDof = (nElement + 1) * 2;
    cout << "number of element = " << nElement << endl;
    cout << "number of dof = " << nDof << endl;

    std::shared_ptr<siconos::algebra::SiconosMatrix> SparseMass(
        new SiconosMatrix(nDof, nDof, siconos::SPARSE, 6 * nDof));
    std::shared_ptr<siconos::algebra::SiconosMatrix> SparseStiffness(
        new SiconosMatrix(nDof, nDof, siconos::SPARSE, 6 * nDof));
    bool lumpedMass = false;
    for (int e = 0; e < nElement; e++) {
      addElementaryStiffnessMatrix(SparseStiffness, e, nDof, l);
      addElementaryMassMatrix(SparseMass, e, nDof, l, lumpedMass);
    }
    // std::cout << "E " << E << std::endl;

    *SparseMass *= rho * S * l / 420.;
    *SparseStiffness *= E * I / (l * l * l);
    cout << "(*SparseMass)(0,0) = " << SparseMass->getValue(0, 0) << endl;
    cout << "(*SparseMass)(0,1) = " << SparseMass->getValue(0, 1) << endl;
    cout << "(*SparseMass)(1,0) = " << SparseMass->getValue(1, 0) << endl;
    cout << "(*SparseMass)(1,1) = " << SparseMass->getValue(1, 1) << endl;

    // SparseMass->display();
    // SparseStiffness->display();
    std::cout << " SparseMass nnz :" << SparseMass->nnz() << std::endl;
    std::cout << " SparseStiffness nnz :" << SparseStiffness->nnz() << std::endl;

    // -- Initial positions and velocities --
    std::shared_ptr<siconos::algebra::SiconosVector> q0(new SiconosVector(nDof, 0.));
    std::shared_ptr<siconos::algebra::SiconosVector> v0(new SiconosVector(nDof, 0.));
    // v0->setValue(0,-.1);
    //  -- The dynamical system --
    auto beam(new LagrangianLinearTIDS(q0, v0, SparseMass));

    // -- Set stiffness matrix (weight) --
    beam->setStiffnessMatrix(SparseStiffness);

    // -- Set external forces (weight) --
    // std::shared_ptr<siconos::algebra::SiconosVector> weight(new SiconosVector(nDof,0.));
    // for (int i =0; i < nDof; i++)
    // {
    //   weight->setValue(i,-g*rho*S/l);
    //   i++;
    // }
    // weight->display();
    Vector weight{nDof};
    weight.setZero();
    beam->etConstantFExt(weight);

    auto bdindex(new IndexInt(2));
    (*bdindex)[0] = nDof - 1;
    (*bdindex)[1] = nDof - 2;

    std::shared_ptr<siconos::algebra::SiconosVector> bdPrescribedVelocity(
        new SiconosVector(2));
    bdPrescribedVelocity->setValue(0, 0.0);
    bdPrescribedVelocity->setValue(1, 0.0);
    auto bd(new BoundaryCondition(bdindex, bdPrescribedVelocity));

    beam->setBoundaryConditions(bd);

    //  Impacting ball

    std::shared_ptr<siconos::algebra::SiconosMatrix> ballMass(
        new SiconosMatrix(1, 1, siconos::SPARSE, 1));
    double ball_mass = 1.0;
    ballMass->setValue(0, 0, ball_mass);

    // std::shared_ptr<siconos::algebra::SiconosMatrix> ballMass(new SiconosMatrix(1,1));
    // ballMass->setValue(0, 0, 1.);

    std::shared_ptr<siconos::algebra::SiconosVector> q0_ball(
        new SiconosVector(1, position_init));
    std::shared_ptr<siconos::algebra::SiconosVector> v0_ball(
        new SiconosVector(1, velocity_init));

    // -- The dynamical system --
    auto ball(new LagrangianLinearTIDS(q0_ball, v0_ball, ballMass));

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.0;

    // Interaction beam-ball
    //
    std::shared_ptr<siconos::algebra::SiconosMatrix> H(new SiconosMatrix(1, nDof + 1));
    (*H)(0, 0) = -1.0;
    (*H)(0, nDof) = 1.0;

    auto nslaw(new NewtonImpactNSL(e));
    auto relation(new LagrangianLinearTIR(H));

    auto inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------
    auto impactingBeam(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    impactingBeam->insertDynamicalSystem(beam);
    impactingBeam->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    impactingBeam->link(inter, beam, ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI(new MoreauJeanOSI(theta, 0.0));
    // OSI->setIsWSymmetricDefinitePositive(true);
    // OSI->setConstraintActivationThreshold(1e-05);

    // -- (2) Time discretisation --
    auto t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    auto osnspb(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    auto s(new TimeStepping(impactingBeam, t, OSI, osnspb));

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    int N = floor((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 17;
    SiconosMatrix dataPlot(N, outputSize);

    std::shared_ptr<siconos::algebra::SiconosVector> q = beam->q();
    std::shared_ptr<siconos::algebra::SiconosVector> v = beam->velocity();
    std::shared_ptr<siconos::algebra::SiconosVector> p = beam->p(1);
    std::shared_ptr<siconos::algebra::SiconosVector> lambda = inter->lambda(1);
    std::shared_ptr<siconos::algebra::SiconosVector> u = inter->y(1);

    int k = 0;

    dataPlot(k, 0) = impactingBeam->t0();

    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*v)(0);
    dataPlot(k, 3) = (*p)(0);
    dataPlot(k, 4) = (*lambda)(0);
    dataPlot(k, 7) = (*q)(nDof - 2);
    dataPlot(k, 8) = (*v)(nDof - 2);

    dataPlot(k, 9) = (*q)((nDof) / 2);
    dataPlot(k, 10) = (*v)((nDof) / 2);

    std::shared_ptr<siconos::algebra::SiconosVector> tmp(new SiconosVector(nDof));

    *tmp = *SparseStiffness * *q;
    double potentialEnergy = 0.5 * q->dot(tmp);
    *tmp = *SparseMass * *v;
    double kineticEnergy = 0.5 * v->dot(tmp);
    double impactEnergy = 0.0;

    dataPlot(k, 5) = potentialEnergy;
    dataPlot(k, 6) = kineticEnergy;
    dataPlot(k, 11) = impactEnergy;

    std::shared_ptr<siconos::algebra::SiconosVector> qBall = ball->q();
    std::shared_ptr<siconos::algebra::SiconosVector> vBall = ball->velocity();
    std::shared_ptr<siconos::algebra::SiconosVector> pBall = ball->p(1);

    dataPlot(k, 12) = (*qBall)(0);
    dataPlot(k, 13) = (*vBall)(0);

    double ballKineticEnergy = 0.5 * (*vBall)(0) * ball_mass * (*vBall)(0);
    dataPlot(k, 16) = ballKineticEnergy;

    //    std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
    //    std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    k++;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    double max_impulse = 0.0;
    //    while (s->nextTime() < T)
    while (k < N) {
      if (k % 100 == 0) std::cout << "k :" << k << "/" << N << std::endl;
      s->computeOneStep();
      // osnspb->display();
      //  std::cout << "\nposition "  << (*q)(0) <<  std::endl;
      //  std::cout << "velocity "  << (*v)(0) <<  std::endl;

      // std::cout << "ball position "  << (*qBall)(0) <<  std::endl;
      // std::cout << "ball velocity "  << (*vBall)(0) <<  std::endl;

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 7) = (*q)(nDof - 2);
      dataPlot(k, 8) = (*v)(nDof - 2);

      dataPlot(k, 9) = (*q)((nDof) / 2);
      dataPlot(k, 10) = (*v)((nDof) / 2);

      *tmp = *SparseStiffness * *q;
      potentialEnergy = 0.5 * q->dot(tmp);
      *tmp = *SparseMass * *v;
      kineticEnergy = 0.5 * v->dot(tmp);

      dataPlot(k, 5) = potentialEnergy;
      dataPlot(k, 6) = kineticEnergy;

      dataPlot(k, 12) = (*qBall)(0);
      dataPlot(k, 13) = (*vBall)(0);
      dataPlot(k, 14) = (*u)(0);

      ballKineticEnergy = 0.5 * (*vBall)(0) * ball_mass * (*vBall)(0);
      // std::cout << "ballKineticEnergy " << ballKineticEnergy << std::endl;
      dataPlot(k, 16) = ballKineticEnergy;

      double u_p_theta = theta * ((*vBall)(0) - (*v)(0)) +
                         (1 - theta) * (dataPlot(k - 1, 13) - dataPlot(k - 1, 2));
      impactEnergy = (u_p_theta) * (*lambda)(0);
      dataPlot(k, 11) = impactEnergy;
      dataPlot(k, 15) = u_p_theta;
      // std::cout <<"\npotentialEnergy ="<<potentialEnergy << std::endl;
      // std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;
      // std::cout <<"impactEnergy ="<<impactEnergy << std::endl;
      // if ((*lambda)(0) >0 )
      //  {
      //   std::cout <<"(*lambda)(0) ="<<(*lambda)(0) << std::endl;
      //   std::cout <<"u_p_theta ="<<u_p_theta << std::endl;
      //   std::cout <<"v beam ="<<(*v)(0) << std::endl;
      //   std::cout <<"v ball ="<<(*vBall)(0) << std::endl;

      //   std::cout <<"v beam before ="<<dataPlot(k-1,2) << std::endl;
      //   std::cout <<"v ball before="<<dataPlot(k-1,13) << std::endl;

      //   getchar();
      // osnspb->display();
      // getchar();
      // }
      if ((*lambda)(0) > 0) {
        max_impulse = fmax(max_impulse, (*lambda)(0));
        // osnspb->display();
        // break;
      }

      s->nextStep();

      k++;
    }

    double totalImpactEnergy = 0.0;
    for (int p = 0; p < k; p++) {
      totalImpactEnergy = totalImpactEnergy + dataPlot(p, 11);
    }
    printf("totalImpactEnergy = %e\n", totalImpactEnergy);
    printf("max_impulse = %e\n", max_impulse);

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << endl;
    ;
    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Computation time : " << elapsed << " ms" << endl;
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("ImpactingEBBeam.dat", "ascii", dataPlot, "noDim");
    // cout << " Comparison with a reference file" << endl;
    // SiconosMatrix dataPlotRef(dataPlot);
    // dataPlotRef.setZero();
    // ioMatrix::read("ImpactingBeam.ref", "ascii", dataPlotRef);

    // double error = (dataPlot - dataPlotRef).normInf() ;
    // cout << "Error = " << error << endl;
    // if(error > 1e-11)
    // {
    //   std::cout << "Warning. The result is rather different from the reference file." <<
    //   std::endl; std::cout << "Error = "<< (dataPlot - dataPlotRef).normInf()<<std::endl;
    //   return 1;
    // }

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
