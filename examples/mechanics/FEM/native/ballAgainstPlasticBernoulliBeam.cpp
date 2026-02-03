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

#include <SiconosFEM.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <cmath>  // for fabs

#include "FENode.hpp"
#include "FiniteElement.hpp"
#include "FiniteElementLinearTIDS.hpp"
#include "FiniteElementModel.hpp"
#include "Material.hpp"
#include "MeshUtils.hpp"
#include "SolverOptions.h"
#include "StressLinearTIR.hpp"

using namespace std;
using namespace siconos::mechanics::fem;
using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;

    auto pStart = std::make_shared<MeshVertex>(0, 0.0, 0.0, 0.0);
    auto pEnd = std::make_shared<MeshVertex>(0, 0.0, 4.0, 0.0);
    // auto pEnd = std::make_shared<MeshVertex>(0, 4.0, 0.0, 0.0);
    // int nbBeams = 1;
    int nbBeams = 4;
    // auto mesh = createBeamMesh(pStart, pEnd, nbBeams, 2);
    int dim = 3;
    int dimBeamDOFs = dim == 2 ? 3 : 6;
    auto mesh = createBeamMesh(pStart, pEnd, nbBeams, dim);
    mesh->display(false);

    int bulk_material_tag = 1;
    int boundary_condition_tag = 2;
    int applied_force_tag = 3;

    double density = 1080.;
    auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(density, 500e7, 1 / 3, 2);
    // auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(
    //     density, 500e5, 1 / 3, 2);

    std::map<unsigned int, std::shared_ptr<siconos::mechanics::fem::Material>> materials = {
        {bulk_material_tag, mat1}, {boundary_condition_tag, mat1}, {applied_force_tag, mat1}};
    auto start = std::chrono::system_clock::now();

    auto beam = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(
        mesh, materials, siconos::algebra::UblasType::SPARSE);
    std::cout << "beam->n(): " << beam->n() << std::endl;
    // int nDof_Beam = beam->n()/2;
    int nDof_Beam = beam->velocityDimension();
    std::cout << "beamTIDS created! ndofs: " << nDof_Beam << std::endl;
    auto femodel = beam->FEModel();
    /*------------------------------------------------- Applied forces  */
    //    (*nodal_forces)(0) = 1e6;
    auto nodal_forces = std::make_shared<Vector>(3);
    nodal_forces->zero();
    (*nodal_forces)(0) = -2e9;
    // beam->applyNodalForces(applied_force_tag, nodal_forces);
    // std::cout << "forces applied!" << std::endl;
    /*------------------------------------------------- Boundary Conditions  */
    /* This part should be hidden in a new BC function for a node number
     * and a dof index. */

    auto node_dof_index = std::make_shared<std::vector<int>>(0);
    node_dof_index->push_back(0);
    node_dof_index->push_back(1);
    node_dof_index->push_back(2);
    beam->applyDirichletBoundaryConditions(boundary_condition_tag, node_dof_index);
    beam->boundaryConditions()->display();
    beam->mass();

    // -------------
    // --- Model ---
    // -------------
    double t0 = 0;               // initial computation time
    double T = 4e-02;            // final computation time
    double h = 1e-05;            // time step
    double theta = 0.5;          // theta for MoreauJeanOSI integrator
    int N = ceil((T - t0) / h);  // Number of time steps

    // unsigned int nDof = 3;  // degrees of freedom for the ball
    unsigned int nslawdimension = 3;  // degrees of freedom for the NS law
    unsigned int nDofBlock = 3;

    double position_init = 0.5;      // initial position for the block.
    double velocity_init = -1000.0;  // initial velocity for the block.

    double R = 3;  // Ball radius
    // double m1 = 300;               // Ball mass
    double m1 = 3500;  // Ball mass

    // auto Mass = std::make_shared<Matrix>(nDof, nDof);

    auto Mass = std::make_shared<Matrix>(nDofBlock, nDofBlock);
    (*Mass)(0, 0) = m1;
    (*Mass)(1, 1) = m1;
    (*Mass)(2, 2) = 2. / 5 * m1 * R * R;

    // -- Initial positions and velocities --
    auto q0 = std::make_shared<Vector>(nDofBlock);
    auto v0 = std::make_shared<Vector>(nDofBlock);

    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;
    double blockHeight = 4.0;
    int Hindex = (int)(blockHeight / (4 / nbBeams)) * dimBeamDOFs;

    // -- The block dynamical system --
    auto block = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, v0, Mass);

    auto beamNSDS = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    beamNSDS->insertDynamicalSystem(beam);
    beamNSDS->insertDynamicalSystem(block);

    // ---------------------------
    // --- Contact Interaction ---
    // ---------------------------

    auto nslawContact =
        std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(0, 0, 0, nslawdimension);
    double conditionning = 1.0;

    auto H_contact = std::make_shared<Matrix>(nslawdimension, nDofBlock + nDof_Beam);
    (*H_contact)(0, Hindex) = -1.0 / conditionning;
    (*H_contact)(1, Hindex + 1) = -1.0 / conditionning;
    if (nslawdimension == 3) (*H_contact)(2, Hindex + 2) = -1.0 / conditionning;
    (*H_contact)(0, (nbBeams + 1) * dimBeamDOFs) = 1.0 / conditionning;
    (*H_contact)(1, (nbBeams + 1) * dimBeamDOFs + 1) = 1.0 / conditionning;
    if (nslawdimension == 3)
      (*H_contact)(2, (nbBeams + 1) * dimBeamDOFs + 2) = 1.0 / conditionning;

    auto b_contact = std::make_shared<Vector>(nslawdimension);

    (*b_contact)(0) = 0.0;
    (*b_contact)(1) = 0.0;
    if (nslawdimension == 3) (*b_contact)(2) = 0.0;

    auto relation_contact =
        std::make_shared<siconos::modeling::LagrangianLinearTIR>(H_contact, b_contact);
    std::cout << "LagrangianLinearTIR created " << std::endl;

    auto inter_contact =
        std::make_shared<siconos::modeling::Interaction>(nslawContact, relation_contact);
    std::cout << "inter_contact created " << std::endl;

    beamNSDS->link(inter_contact, beam, block);

    // ------------------------------
    // --- Plasticity Interaction ---
    // ------------------------------

    // double c = 1.0e8;  // Used only to specify elastic bounds
    double c = 1.0e9;  // Used only to specify elastic bounds

    double phi = M_PI / 4;  // Not used at the moment

    auto nslawPlasticity =
        std::make_shared<siconos::modeling::MohrCoulombPlasticityNSL>(c, phi, nslawdimension);

    auto initialStress_gap = std::make_shared<Vector>(nslawdimension, 0);

    for (int elNum = 3; elNum < 4; elNum++) {
      auto H = std::make_shared<Matrix>(nslawdimension, beam->stressDimension());
      auto Hminus = std::make_shared<Matrix>(nslawdimension, beam->stressDimension());
      auto H2 = std::make_shared<Matrix>(nslawdimension, beam->stressDimension());
      auto Hminus2 = std::make_shared<Matrix>(nslawdimension, beam->stressDimension());

      (*H)(0, elNum * 3 + 1) = 1.0;
      (*H)(1, elNum * 3 + 2) = 1.0;
      if (nslawdimension == 3) (*H)(2, elNum * 3) = 1.0;
      (*Hminus)(0, elNum * 3 + 1) = -1.0;
      (*Hminus)(1, elNum * 3 + 2) = -1.0;
      if (nslawdimension == 3) (*Hminus)(2, elNum * 3) = -1.0;
      (*H2)(0, elNum * 3 + 2) = 1.0;
      (*H2)(1, elNum * 3) = 1.0;
      if (nslawdimension == 3) (*H2)(2, elNum * 3 + 1) = 1.0;
      (*Hminus2)(0, elNum * 3 + 2) = -1.0;
      (*Hminus2)(1, elNum * 3) = -1.0;
      if (nslawdimension == 3) (*Hminus2)(2, elNum * 3 + 1) = -1.0;

      auto relation =
          std::make_shared<siconos::modeling::StressLinearTIR>(H, initialStress_gap);
      auto inter = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relation);
      auto relationMinus =
          std::make_shared<siconos::modeling::StressLinearTIR>(Hminus, initialStress_gap);
      auto interMinus =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationMinus);
      auto relation2 =
          std::make_shared<siconos::modeling::StressLinearTIR>(H2, initialStress_gap);
      auto inter2 =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relation2);
      auto relationMinus2 =
          std::make_shared<siconos::modeling::StressLinearTIR>(Hminus2, initialStress_gap);
      auto interMinus2 =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationMinus2);

      // beamNSDS->link(inter, beam);
      // beamNSDS->link(interMinus, beam);
      // beamNSDS->link(inter2, beam);
      // beamNSDS->link(interMinus2, beam);
    }
    std::cout << "Out loop plsticity" << std::endl;

    // ---------------------------
    // --- Integration and solver Setup ---
    // ---------------------------

    auto OSI = std::make_shared<siconos::integrators::MoreauJeanGOSI>(theta);
    // OSI->setIsWSymmetricDefinitePositive(true);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // auto osnspb =
    // std::make_shared<siconos::nonsmooth_formulations::GlobalFrictionContact>(nslawdimension
    // , SICONOS_FRICTION_2D_NSGS);
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::GlobalFrictionContact>(
        nslawdimension, SICONOS_GLOBAL_FRICTION_3D_NSGS_WR);
    osnspb->setNumericsVerboseLevel(1);

    osnspb->numericsSolverOptions()->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] =
        SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL;
    // osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1.0e-1;

    // auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    std::cout << "beam->stressDimension() - 2: " << beam->stressDimension() - 2 << std::endl;
    auto s = std::make_shared<siconos::simulation::TimeStepping>(beamNSDS, t, OSI, osnspb);

    std::string SOFAfilename = "beam";
    int k = 1;
    start = std::chrono::system_clock::now();

    auto qball = block->q();
    auto vball = block->velocity();
    auto q = beam->q();
    auto v = beam->velocity();
    auto sigma = beam->stress();
    auto p1 = block->p(1);
    auto epsilonPointP = beam->plasticRate_read();
    auto epsilonPointPOld =
        std::make_shared<siconos::algebra::SiconosVector>(epsilonPointP.size());
    epsilonPointPOld->fill(0);
    auto diff = std::make_shared<siconos::algebra::SiconosVector>(epsilonPointP.size());
    auto epsilonE = std::make_shared<siconos::algebra::SiconosVector>(epsilonPointP.size());
    auto epsilon = std::make_shared<siconos::algebra::SiconosVector>(epsilonPointP.size());
    auto sigmaOld = std::make_shared<siconos::algebra::SiconosVector>(beam->stressDimension());
    sigmaOld->fill(0);
    auto vMass = std::make_shared<siconos::algebra::SiconosVector>(beam->velocity()->size());
    auto vBallMass =
        std::make_shared<siconos::algebra::SiconosVector>(block->velocity()->size());
    auto vTotal = std::make_shared<siconos::algebra::SiconosVector>(beam->velocityDimension() +
                                                                    block->velocity()->size());
    auto vContact = std::make_shared<siconos::algebra::SiconosVector>(nslawdimension);
    auto vContactOld = std::make_shared<siconos::algebra::SiconosVector>(nslawdimension);

    vTotal->segment(0, beam->dimension()) = *v;
    vTotal->segment(beam->dimension(), block->velocity()->size()) = *vball;
    // vTotal = [vTree; vBlock]

    siconos::algebra::prod(*H_contact, *vTotal, *vContactOld, true);  // Init vContactOld

    prepareWriteBeamPositionforSOFA(SOFAfilename + ".state");
    writeBeamPositionforSOFA(mesh, femodel, q, SOFAfilename + ".state", 0);

    prepareWriteBlockPositionforSOFA(SOFAfilename + "_block.state");

    (*q0)(1) = blockHeight;
    writeBlockPositionforSOFA(q0, SOFAfilename + "_block.state", 0);

    unsigned int outputSize = 11;
    Matrix dataPlot(N + 1, outputSize);

    dataPlot(0, 0) = beamNSDS->t0();
    dataPlot(0, 1) = (*sigma)(beam->stressDimension() - 2);
    dataPlot(0, 2) = (*sigma)(beam->stressDimension() - 1);
    dataPlot(0, 3) = (*epsilonPointP)(beam->stressDimension() - 2);
    dataPlot(0, 4) = (*epsilonPointP)(beam->stressDimension() - 1);
    dataPlot(0, 5) = 0;
    dataPlot(0, 6) = 0;
    dataPlot(0, 7) = 0;
    dataPlot(0, 8) = 0;
    dataPlot(0, 9) = dataPlot(0, 4);
    dataPlot(0, 10) = (*p1)(0);

    double energieCinetiqueBlock, dissipationContact = 0, dissipationPlastique = 0,
                                  energieMecanique, energieCinetiqueTree;
    auto filename = siconos::mechanics::fem::prepareWriteDisplacementforPython("beam");

    while (s->hasNextEvent()) {
      std::cout << "Step number " << k << std::endl;
      s->computeOneStep();
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*sigma)(beam->stressDimension() - 2);
      dataPlot(k, 2) = (*sigma)(beam->stressDimension() - 1);
      dataPlot(k, 3) = (*epsilonPointP)(beam->stressDimension() - 2);
      siconos::algebra::prod(*(block->mass()), *vball, *vBallMass, true);
      energieCinetiqueBlock = (1.0 / 2.0) * siconos::algebra::inner_prod(*vBallMass, *vball);
      siconos::algebra::prod(*(beam->mass()), *v, *vMass, true);
      energieCinetiqueTree = (1.0 / 2.0) * siconos::algebra::inner_prod(*vMass, *v);
      siconos::algebra::sub(*epsilonPointP, *epsilonPointPOld, *diff);
      epsilonPointPOld->fill(0);
      siconos::algebra::axpy(1.0, *epsilonPointP, *epsilonPointPOld);

      siconos::algebra::axpby(theta, *sigma, (1 - theta), *sigmaOld);

      dissipationPlastique += siconos::algebra::inner_prod(*sigmaOld, *diff);
      sigmaOld->fill(0);
      siconos::algebra::axpy(1.0, *sigma, *sigmaOld);
      siconos::algebra::prod(*(beam->B()), *q, *epsilonE, true);
      siconos::algebra::sub(*epsilonE, *epsilonPointP, *epsilon);
      energieMecanique = (1.0 / 2.0) * siconos::algebra::inner_prod(*sigma, *epsilon);

      vTotal->segment(0, beam->dimension()) = *v;
      vTotal->segment(beam->dimension(), block->velocity()->size()) = *vball;
      // vTotal = [vTree; vBlock]

      siconos::algebra::prod(*H_contact, *vTotal, *vContact, true);
      siconos::algebra::axpby(theta, *vContact, (1 - theta), *vContactOld);
      dissipationContact -= siconos::algebra::inner_prod(*vContactOld, *p1);
      siconos::algebra::scal(1.0, *vContact, *vContactOld, true);  // reset vContactOld

      dataPlot(k, 4) = (*epsilonPointP)(beam->stressDimension() - 1);
      dataPlot(k, 5) = energieCinetiqueTree;
      dataPlot(k, 6) = energieMecanique;
      dataPlot(k, 7) = dissipationPlastique;
      dataPlot(k, 8) = dissipationContact;
      dataPlot(k, 9) = energieCinetiqueBlock + energieMecanique + energieCinetiqueTree +
                       dissipationPlastique + dissipationContact;
      dataPlot(k, 10) = (*p1)(0);

      writeBeamPositionforSOFA(mesh, femodel, q, SOFAfilename + ".state", k * 0.01);
      writeBlockPositionforSOFA(qball, SOFAfilename + "_block.state", k * 0.01);

      if (k % 1 == 0) writeDisplacementforPython(mesh, femodel, q, filename);

      s->nextStep();
      k++;
      siconos::tools::progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    siconos::algebra::io::write("ballAgainstPlasticBernoulliBeam.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    std::cout << "Done! " << std::endl;

    double error = 0.0, eps = 1e-12;
    // if ((error = siconos::algebra::io::compareRefFile(
    //          dataPlot, "bernoulliBeam.ref", eps)) >= eps)
    //   return 1;

    return 0;

  } catch (...) {
    std::cerr << "Exception caught in ballAgainstBernoulliBeam.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
