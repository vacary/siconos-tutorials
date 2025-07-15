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
#include "FiniteElementModel.hpp"
#include "MeshUtils.hpp"
#include "FiniteElementLinearTIDS.hpp"
#include "Material.hpp"
#include "StressLinearTIR.hpp"

using namespace std;
using namespace siconos::mechanics::fem;
using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;

    auto p0 = std::make_shared<MVertex>(0, 0.0, 0.0, 0.0);
    auto pEnd = std::make_shared<MVertex>(0, 0.0, 4.0, 0.0);
    int nbBeams = 4;
    int dim = 2;
    auto mesh = createBeamMesh(p0, pEnd, nbBeams, 2);
    mesh->display(false);

            // siconos::mechanics::fem::writeMeshforPython(mesh);

    int bulk_material_tag = 1;
    int boundary_condition_tag = 2;
    int applied_force_tag = 3;

            // std::shared_ptr<Material> mat1 = std::make_shared<Material>(1, 8*36/5.,
            // 1/5.); // material for  triangle_felippa.msh
    double density = 1080.;
    auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(
        density, 500e7, 1 / 3);
    // auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(
    //     1.0, 1.0, 1 / 3);
    std::map<unsigned int, std::shared_ptr<siconos::mechanics::fem::Material>>
        materials = {{bulk_material_tag, mat1},{boundary_condition_tag,mat1},{applied_force_tag,mat1}};
    auto start = std::chrono::system_clock::now();
    std::cout << "Creating beamTIDS... " << std::endl;
    // auto beam = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(
    //     mesh, materials, siconos::algebra::UblasType::SPARSE);
    // int nDof_Beam = beam->velocityDimension();

    // auto beam = std::make_shared<siconos::mechanics::fem::FiniteElementLinearTIDS>(
    //     mesh, materials, siconos::algebra::UblasType::SPARSE);
    auto beam = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(
        mesh, materials, siconos::algebra::UblasType::SPARSE);
    std::cout << "beam->n(): " << beam->n() << std::endl;
    // int nDof_Beam = beam->n()/2;
    int nDof_Beam = beam->velocityDimension();
    std::cout << "beamTIDS created! ndofs: " << nDof_Beam << std::endl;
    auto femodel = beam->FEModel();
    /*------------------------------------------------- Applied forces  */
                                                                            //    (*nodal_forces)(0) = 1e6;
    std::cout << "Applying forces :" << std::endl;

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
    std::cout << "Applying BCs :" << std::endl;
    beam->applyDirichletBoundaryConditions(boundary_condition_tag,
                                           node_dof_index);
    std::cout << "BCs :" << std::endl;
    beam->boundaryConditions()->display();

            // -------------
            // --- Model ---
            // -------------
    double t0 = 0;       // initial computation time
    double T = 4e-02;    // final computation time
    double h = 1e-05;    // time step
    double theta = 1.0;  // theta for MoreauJeanOSI integrator
    int N = ceil((T - t0) / h);  // Number of time steps



    unsigned int nDof = 3;  // degrees of freedom for the ball
    double position_init = 0.5;  // initial position for the block.
    double velocity_init = -1000.0;  // initial velocity for the block.

    double R = 3;              // Ball radius
    double m1 = 300;               // Ball mass
    auto Mass = std::make_shared<Matrix>(nDof, nDof);
    (*Mass)(0, 0) = m1;
    (*Mass)(1, 1) = m1;
    (*Mass)(2, 2) = 2. / 5 * m1 * R * R;

            // -- Initial positions and velocities --
    auto q0 = std::make_shared<Vector>(nDof);
    auto v0 = std::make_shared<Vector>(nDof);
    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;

            // -- The block dynamical system --
    auto block = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, v0, Mass);




    auto beamNSDS =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

            // add the dynamical system in the non smooth dynamical system
    beamNSDS->insertDynamicalSystem(beam);
    beamNSDS->insertDynamicalSystem(block);




    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(0, 0, 0, 3);
    // double e = 0;
    // auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    // auto H_bb = std::make_shared<Matrix>(1, nDof + nDof_Beam);
    // (*H_bb)(0, 9) = -1.0;
    // // (*H_bb)(1, 10) = -1.0;
    // // (*H_bb)(2, 11) = -1.0;
    // (*H_bb)(0, 15) = 1.0;
    // // (*H_bb)(1, 16) = 1.0;
    // // (*H_bb)(2, 17) = 1.0;

            // auto b_bb = std::make_shared<Vector>(1);
            // (*b_bb)(0) = 0.0;
            // // (*b_bb)(1) = 0.0;
            // // (*b_bb)(2) = 0.0;

    // double conditionning = 2.500000e+04;
    double conditionning = 1.0;
    auto H_bb = std::make_shared<Matrix>(nDof, nDof + nDof_Beam);
    (*H_bb)(0, 12) = -1.0/conditionning;
    (*H_bb)(1, 13) = -1.0/conditionning;
    (*H_bb)(2, 14) = -1.0/conditionning;
    (*H_bb)(0, 15) = 1.0/conditionning;
    (*H_bb)(1, 16) = 1.0/conditionning;
    (*H_bb)(2, 17) = 1.0/conditionning;

    auto b_bb = std::make_shared<Vector>(3);
    (*b_bb)(0) = 0.0;
    (*b_bb)(1) = 0.0;
    (*b_bb)(2) = 0.0;

    auto relation_bb = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H_bb, b_bb);

    auto inter_bb = std::make_shared<siconos::modeling::Interaction>(nslaw, relation_bb);


    beamNSDS->link(inter_bb, beam, block);


    double c = 1.0e8;
    double phi = M_PI/4;
    int nslawdimension = 3;

    auto nslawPlasticity = std::make_shared<siconos::modeling::MohrCoulombPlasticityNSL>(c,phi,nslawdimension);

    auto initialStress_gap = std::make_shared<Vector>(nslawdimension, 0);
    int elcount = 2;
    std::cout << "Stress dimension: " << beam->stressDimension() << std::endl;
    for (auto& el : femodel->elements()) {
      auto Htension = std::make_shared<Matrix>(nslawdimension, beam->stressDimension());
      auto HtensionMinus = std::make_shared<Matrix>(nslawdimension, beam->stressDimension());
      auto H = std::make_shared<Matrix>(nslawdimension, beam->stressDimension());
      auto Hminus = std::make_shared<Matrix>(nslawdimension, beam->stressDimension());
      auto H2 = std::make_shared<Matrix>(nslawdimension, beam->stressDimension());
      auto Hminus2 = std::make_shared<Matrix>(nslawdimension, beam->stressDimension());

      el->display();
      std::cout << "elcount " << elcount << std::endl;
      (*Htension)(0, elcount*3) = 1.0;
      (*Htension)(1, elcount*3+1) = 1.0;
      if (nslawdimension == 3)
        (*H)(2, elcount*3+2) = 1.0;
      (*HtensionMinus)(0, elcount*3) = -1.0;
      (*HtensionMinus)(1, elcount*3+1) = -1.0;
      if (nslawdimension == 3)
        (*HtensionMinus)(2, elcount*3+1) = -1.0;

      (*H)(0, elcount*3+1) = 1.0;
      (*H)(1, elcount*3+2) = 1.0;
      if (nslawdimension == 3)
        (*H)(2, elcount*3+2) = 1.0;
      (*Hminus)(0, elcount*3+1) = -1.0;
      (*Hminus)(1, elcount*3+2) = -1.0;
      if (nslawdimension == 3)
        (*Hminus)(2, elcount*3+1) = -1.0;
      (*H2)(0, elcount*3+2) = 1.0;
      (*H2)(1, elcount*3) = 1.0;
      if (nslawdimension == 3)
        (*H2)(2, elcount*3+2) = 1.0;
      (*Hminus2)(0, elcount*3+2) = -1.0;
      (*Hminus2)(1, elcount*3) = -1.0;
      if (nslawdimension == 3)
        (*Hminus2)(2, elcount*3+1) = -1.0;

      // auto relationTension = std::make_shared<siconos::modeling::StressLinearTIR>(Htension, initialStress_gap);
      // auto interTension = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationTension);
      // auto relationTensionMinus = std::make_shared<siconos::modeling::StressLinearTIR>(HtensionMinus, initialStress_gap);
      // auto interTensionMinus = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationTensionMinus);

      auto relation = std::make_shared<siconos::modeling::StressLinearTIR>(H, initialStress_gap);
      auto inter = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relation);
      auto relationMinus = std::make_shared<siconos::modeling::StressLinearTIR>(Hminus, initialStress_gap);
      auto interMinus = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationMinus);
      auto relation2 = std::make_shared<siconos::modeling::StressLinearTIR>(H2, initialStress_gap);
      auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relation2);
      auto relationMinus2 = std::make_shared<siconos::modeling::StressLinearTIR>(Hminus2, initialStress_gap);
      auto interMinus2 = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationMinus2);

      // beamNSDS->link(interTension, beam);
      // beamNSDS->link(interTensionMinus, beam);
      beamNSDS->link(inter, beam);
      beamNSDS->link(interMinus, beam);
      beamNSDS->link(inter2, beam);
      beamNSDS->link(interMinus2, beam);

      elcount++;
      if (elcount == 3)
        break;
    }



    auto OSI = std::make_shared<siconos::integrators::MoreauJeanGOSI>(theta);
    // auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
    OSI->setIsWSymmetricDefinitePositive(true);

            // -- (2) Time discretisation --

    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // auto osnspb = std::make_shared<siconos::nonsmooth_formulations::GlobalFrictionContact>(nslawdimension , SICONOS_FRICTION_2D_NSGS);
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::GlobalFrictionContact>(nslawdimension,SICONOS_GLOBAL_FRICTION_3D_NSGS_WR);

    // auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();


            // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(beamNSDS, t,OSI,osnspb);
    // s->insertIntegrator(OSI);

    std::string SOFAfilename = "beam.state";
    int k = 1;
    start = std::chrono::system_clock::now();
    unsigned int outputSize = 9;
    Matrix dataPlot(N + 1, outputSize);
    auto qball = block->q();
    auto q = beam->q();
    auto v = beam->velocity();
    auto sigma = beam->stress();
    auto p1 = block->p(1);

    prepareWriteBeamPositionforSOFA(SOFAfilename);
    writeBeamPositionforSOFA(mesh, femodel, q, SOFAfilename,0);

    Matrix pos(N + 1, 10);
    int vcnt = 0;
    // for (auto v : mesh->vertices()) {
    //   std::cout << v->num() << std::endl;
    //   pos(0,(vcnt)*2) = v->x() + (*q)((vcnt)*3);
    //   pos(0,(vcnt)*2+1) = v->y() + (*q)((vcnt)*3+1);
    //   vcnt++;
    // }

    dataPlot(0, 0) = beamNSDS->t0();
    dataPlot(0, 1) = (*q)(beam->dimension() - 6);
    dataPlot(0, 2) = (*q)(beam->dimension() - 5);
    dataPlot(0, 3) = (*q)(beam->dimension() - 4);
    dataPlot(0, 4) = (*sigma)(beam->stressDimension() - 6);
    dataPlot(0, 5) = (*sigma)(beam->stressDimension() - 5);
    dataPlot(0, 6) = (*sigma)(beam->stressDimension() - 4);
    dataPlot(0, 7) = (*qball)(0);
    dataPlot(0, 8) = (*p1)(0);
    auto filename =
        siconos::mechanics::fem::prepareWriteDisplacementforPython("beam");
    auto tensorFilename =
        siconos::mechanics::fem::prepareWriteTensorforPython("beamTensor","tensor");

    std::cout << "About to start simulation" << std::endl;
    while (s->hasNextEvent()) {
      std::cout << "Step number " << k << std::endl;
      s->computeOneStep();
      std::cout << "computed! " << std::endl;
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(beam->dimension() - 6);
      dataPlot(k, 2) = (*q)(beam->dimension() - 5);
      dataPlot(k, 3) = (*q)(beam->dimension() - 4);
      dataPlot(k, 4) = (*sigma)(beam->stressDimension() - 6);
      dataPlot(k, 5) = (*sigma)(beam->stressDimension() - 5);
      dataPlot(k, 6) = (*sigma)(beam->stressDimension() - 4);
      dataPlot(k, 7) = (*qball)(0);
      dataPlot(k, 8) = (*p1)(0);
      vcnt = 0;
      // for (auto v : mesh->vertices()) {
      //   std::cout << v->num() << std::endl;

              //   pos(k,(vcnt)*2) = v->x() + (*q)((vcnt)*3);
              //   pos(k,(vcnt)*2+1) = v->y() + (*q)((vcnt)*3+1);
              //   vcnt++;
              // }
      std::cout << "writing beam pos for SOFA... " << std::endl;
      q->display();
      qball->display();
      sigma->display();
      writeBeamPositionforSOFA(mesh, femodel, q, SOFAfilename,k*0.01);
      if (k % 1 == 0) writeDisplacementforPython(mesh, femodel, q, filename);
      if (k % 1 == 0) writeTensorforPython(femodel, sigma, tensorFilename, "tensor");

      std::cout << "Done! " << std::endl;
      s->nextStep();
      k++;

      siconos::tools::progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                       .count();


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
