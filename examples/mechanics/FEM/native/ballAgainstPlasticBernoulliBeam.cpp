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

#include <FEM-GlobalFrictionContact.hpp>
#include <FEM-MoreauJeanGOSI.hpp>
#include <SiconosKernel.hpp>
#include <StressLinearTIR.hpp>
#include <chrono>
#include <cmath>  // for fabs

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;

    int nb_elements = 4;
    int dim = 3;
    siconos::algebra::SiconosVector3 coords_start{0., 0., 0.};
    siconos::algebra::SiconosVector3 coords_end{0., 4., 0.};

    int dimBeamDOFs = dim == 2 ? 3 : 6;
    auto mesh =
        siconos::mechanics::fem::createBeamMesh(coords_start, coords_end, nb_elements, dim);
    mesh->display(false);

    siconos::mechanics::fem::Tags tags;
    tags[siconos::mechanics::fem::MeshTags::bulk_material] = 1;
    tags[siconos::mechanics::fem::MeshTags::boundary_conditions] = 2;
    tags[siconos::mechanics::fem::MeshTags::applied_forces] = 3;

    siconos::mechanics::fem::Material mat{1080, 500e7, 1. / 3, 2.};

    // Same material for all tags.
    // Same material for all tags.
    std::map<int, const siconos::mechanics::fem::Material> materials = {
        {1, mat}, {2, mat}, {3, mat}};

    // Create the beam
    auto beam = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(mesh, materials);

    // Apply nodal forces
    if (tags.find(siconos::mechanics::fem::MeshTags::applied_forces) != tags.end()) {
      siconos::algebra::SiconosVector nodal_forces{3};
      nodal_forces << -2e9, 0., 0.;

      beam->applyNodalForces(tags[siconos::mechanics::fem::MeshTags::applied_forces],
                             nodal_forces);
    }
    // Boundary Conditions

    if (tags.find(siconos::mechanics::fem::MeshTags::boundary_conditions) != tags.end()) {
      std::vector<int> bc_dof_index(3);
      bc_dof_index[0] = 0;
      bc_dof_index[1] = 1;
      bc_dof_index[1] = 2;
      beam->applyDirichletBoundaryConditions(
          tags[siconos::mechanics::fem::MeshTags::boundary_conditions], bc_dof_index);
    }

    std::cout << beam->stressDimension() << "\n";

    // ------- Now the block -------

    double R = 3;      // Ball radius
    double m1 = 3500;  // Ball mass
    siconos::algebra::SiconosSparseMatrix mass{3, 3};
    mass.insert(0, 0) = m1;
    mass.insert(1, 1) = m1;
    mass.insert(2, 2) = 2. / 5 * m1 * R * R;
    mass.makeCompressed();

    // -- Initial positions and velocities --
    double position_init = 0.5;      // initial position for the block.
    double velocity_init = -1000.0;  // initial velocity for the block.
    siconos::algebra::SiconosVector3 q0;
    siconos::algebra::SiconosVector3 v0;
    q0 << position_init, 0., 0.;
    v0 << velocity_init, 0., 0.;

    auto block = std::make_shared<siconos::modeling::LagrangianSparseLinearTIDS>(
        q0, v0, mass, siconos::algebra::copy_t);

    // ------- Interactions -------

    auto ndof_beam = beam->dimension();
    auto ndof_ball = block->dimension();
    // -- nslaw --
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(0, 0, 0, 3);
    // Interaction ball-floor
    siconos::algebra::SiconosDenseMatrix H_bb{ndof_ball, ndof_ball + ndof_beam};
    double conditionning = 1.0;

    H_bb.setZero();
    double blockHeight = 4.0;
    size_t Hindex = (size_t)(blockHeight / (4 / nb_elements)) * dimBeamDOFs;

    H_bb(0, Hindex) = -1.0 / conditionning;
    H_bb(1, Hindex + 1) = -1.0 / conditionning;
    H_bb(2, Hindex + 2) = -1.0 / conditionning;
    int pos = (nb_elements + 1) * dimBeamDOFs;
    H_bb(0, pos) = 1.0 / conditionning;
    H_bb(1, pos + 1) = 1.0 / conditionning;
    H_bb(2, pos + 2) = 1.0 / conditionning;

    auto relation_bb = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H_bb);
    auto inter_bb = std::make_shared<siconos::modeling::Interaction>(nslaw, relation_bb);

    // ------- NSDS -------
    double t0 = 0;     // initial computation time
    double T = 4e-02;  // final computation time
    auto beamNSDS = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    beamNSDS->insertDynamicalSystem(beam);
    beamNSDS->insertDynamicalSystem(block);
    beamNSDS->link(inter_bb, beam, block);

    // unsigned int nDof = 3;  // degrees of freedom for the ball
    siconos::algebra::Index nslawdimension = 3;  // degrees of freedom for the NS law

    // ------------------------------
    // --- Plasticity Interaction ---
    // ------------------------------

    // double c = 1.0e8;  // Used only to specify elastic bounds
    double c = 1.0e9;  // Used only to specify elastic bounds

    double phi = M_PI / 4;  // Not used at the moment

    auto nslawPlasticity =
        std::make_shared<siconos::modeling::MohrCoulombPlasticityNSL>(c, phi, nslawdimension);

    std::vector<siconos::algebra::SiconosDenseMatrix> Hv;
    std::vector<siconos::algebra::SiconosDenseMatrix> Hminus_v;
    std::vector<siconos::algebra::SiconosDenseMatrix> H2v;
    std::vector<siconos::algebra::SiconosDenseMatrix> Hminus2_v;
    for (int elNum = 3; elNum < 4; elNum++) {
      Hv.emplace_back(nslawdimension, beam->stressDimension());

      Hminus_v.emplace_back(nslawdimension, beam->stressDimension());
      H2v.emplace_back(nslawdimension, beam->stressDimension());
      Hminus2_v.emplace_back(nslawdimension, beam->stressDimension());

      Hv.back().setZero();
      Hv.back()(0, elNum * 3 + 1) = 1.0;
      Hv.back()(1, elNum * 3 + 2) = 1.0;
      if (nslawdimension == 3) Hv.back()(2, elNum * 3) = 1.0;
      Hminus_v.back()(0, elNum * 3 + 1) = -1.0;
      Hminus_v.back()(1, elNum * 3 + 2) = -1.0;
      if (nslawdimension == 3) Hminus_v.back()(2, elNum * 3) = -1.0;
      H2v.back()(0, elNum * 3 + 2) = 1.0;
      H2v.back()(1, elNum * 3) = 1.0;
      if (nslawdimension == 3) H2v.back()(2, elNum * 3 + 1) = 1.0;
      Hminus2_v.back()(0, elNum * 3 + 2) = -1.0;
      Hminus2_v.back()(1, elNum * 3) = -1.0;
      if (nslawdimension == 3) Hminus2_v.back()(2, elNum * 3 + 1) = -1.0;

      auto relation = std::make_shared<siconos::mechanics::fem::StressLinearTIR>(Hv.back());
      auto inter = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relation);
      auto relationMinus =
          std::make_shared<siconos::mechanics::fem::StressLinearTIR>(Hminus_v.back());
      auto interMinus =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationMinus);
      auto relation2 = std::make_shared<siconos::mechanics::fem::StressLinearTIR>(H2v.back());
      auto inter2 =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relation2);
      auto relationMinus2 =
          std::make_shared<siconos::mechanics::fem::StressLinearTIR>(Hminus2_v.back());
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
    // ------- Simulation -------
    // -- (1) integrator --
    double theta = 0.5;  // theta for MoreauJeanOSI integrator
    auto osi_fem =
        std::make_shared<siconos::mechanics::fem::integrators::MoreauJeanGOSI>(theta);
    // osi_fem->setIsWSymmetricDefinitePositive(true);

    // -- (2) Time discretisation --
    double h = 1e-05;  // time step
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem --
    auto osnspb = std::make_shared<
        siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact>(
        nslawdimension, SICONOS_GLOBAL_FRICTION_3D_NSGS_WR);

    osnspb->setNumericsVerboseLevel(1);

    osnspb->numericsSolverOptions()->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] =
        SICONOS_NSGS_ERROR_EVALUATION_FULL;
    // osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1.0e-1;

    // auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    std::cout << "beam->stressDimension() - 2: " << beam->stressDimension() - 2 << std::endl;
    auto s = std::make_shared<siconos::simulation::TimeStepping>(beamNSDS, t, osi_fem, osnspb);

    int N = ceil((T - t0) / h);  // Number of time steps
    std::string SOFAfilename = "beam";
    int k = 1;

    siconos::algebra::Index outputSize = 11;
    siconos::algebra::SiconosDenseMatrix dataPlot{N + 1, outputSize};

    auto start = std::chrono::system_clock::now();

    auto qball = block->q_read();
    auto q = beam->q_read();
    auto v = beam->velocity_read();
    auto p1 = block->p_read(1);

    auto vball = block->velocity_read();
    auto sigma = beam->stress_read();
    auto epsilonPointP = beam->plasticRate_read();

    siconos::algebra::SiconosVector vTotal{beam->dimension() + block->velocity()->size()};
    vTotal.segment(0, beam->dimension()) = v;
    vTotal.segment(beam->dimension(), block->velocity()->size()) = vball;
    // vTotal = [vTree; vBlock]
    siconos::algebra::SiconosVector vContactOld = H_bb * vTotal;

    auto femodel = beam->FEModel();
    siconos::mechanics::fem::prepareWriteBeamPositionforSOFA(SOFAfilename + ".state");
    siconos::mechanics::fem::writeBeamPositionforSOFA(*mesh, *femodel, q,
                                                      SOFAfilename + ".state", 0);

    siconos::mechanics::fem::prepareWriteBlockPositionforSOFA(SOFAfilename + "_block.state");

    q0(1) = blockHeight;
    siconos::mechanics::fem::writeBlockPositionforSOFA(q0, SOFAfilename + "_block.state", 0);

    dataPlot(0, 0) = beamNSDS->t0();
    dataPlot(0, 1) = sigma(beam->stressDimension() - 2);
    dataPlot(0, 2) = sigma(beam->stressDimension() - 1);
    dataPlot(0, 3) = epsilonPointP(beam->stressDimension() - 2);
    dataPlot(0, 4) = epsilonPointP(beam->stressDimension() - 1);
    dataPlot(0, 5) = 0;
    dataPlot(0, 6) = 0;
    dataPlot(0, 7) = 0;
    dataPlot(0, 8) = 0;
    dataPlot(0, 9) = dataPlot(0, 4);
    dataPlot(0, 10) = p1(0);

    double dissipationContact = 0, dissipationPlastique = 0, energieMecanique,
           energieCinetiqueTree;
    auto filename = siconos::mechanics::fem::prepareWriteDisplacementforPython("beam");
    siconos::algebra::SiconosVector epsilonPointPOld{epsilonPointP.size()};
    epsilonPointPOld.setZero();
    siconos::algebra::SiconosVector sigmaOld{beam->stressDimension()};
    sigmaOld.setZero();

    while (s->hasNextEvent()) {
      std::cout << "Step number " << k << std::endl;
      s->computeOneStep();
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = sigma(beam->stressDimension() - 2);
      dataPlot(k, 2) = sigma(beam->stressDimension() - 1);
      dataPlot(k, 3) = epsilonPointP(beam->stressDimension() - 2);
      auto energieCinetiqueBlock = 0.5 * (block->mass() * vball).dot(vball);

      energieCinetiqueTree = 0.5 * (beam->mass() * v).dot(v);
      auto diff = epsilonPointP - epsilonPointPOld;
      epsilonPointPOld = epsilonPointP;
      sigmaOld = theta * sigma + (1. - theta) * sigmaOld;

      dissipationPlastique += sigmaOld.dot(diff);
      sigmaOld = sigma;
      auto epsilonE = beam->BMatrix() * q;
      auto epsilon = epsilonE - epsilonPointP;
      energieMecanique = 0.5 * sigma.dot(epsilon);

      vTotal.segment(0, beam->dimension()) = v;
      vTotal.segment(beam->dimension(), block->velocity()->size()) = vball;
      // vTotal = [vTree; vBlock]
      auto vContact = H_bb * vTotal;
      vContactOld = theta * vContact + (1. - theta) * vContactOld;
      dissipationContact -= vContactOld.dot(p1);
      vContactOld = vContact;

      dataPlot(k, 4) = epsilonPointP(beam->stressDimension() - 1);
      dataPlot(k, 5) = energieCinetiqueTree;
      dataPlot(k, 6) = energieMecanique;
      dataPlot(k, 7) = dissipationPlastique;
      dataPlot(k, 8) = dissipationContact;
      dataPlot(k, 9) = energieCinetiqueBlock + energieMecanique + energieCinetiqueTree +
                       dissipationPlastique + dissipationContact;
      dataPlot(k, 10) = p1(0);

      siconos::mechanics::fem::writeBeamPositionforSOFA(*mesh, *femodel, q,
                                                        SOFAfilename + ".state", k * 0.01);
      siconos::mechanics::fem::writeBlockPositionforSOFA(qball, SOFAfilename + "_block.state",
                                                         k * 0.01);

      if (k % 1 == 0)
        siconos::mechanics::fem::writeDisplacementforPython(*mesh, *femodel, q, filename);

      s->nextStep();
      k++;

      siconos::tools::progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    siconos::algebra::io::write("ballAgainstPlasticBernoulliBeam.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

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
