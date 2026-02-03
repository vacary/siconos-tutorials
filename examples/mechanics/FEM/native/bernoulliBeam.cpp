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
#include "FiniteElementModel.hpp"
#include "MeshUtils.hpp"
#include "FiniteElementLinearTIDS.hpp"
#include "Material.hpp"

using namespace std;
using namespace siconos::mechanics::fem;
using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;

    auto v0 = std::make_shared<MeshVertex>(0, 0.0, 0.0, 0.0);
    auto vEnd = std::make_shared<MeshVertex>(0, 0.0, 4.0, 0.0);
    int nbBeams = 4;
    int dim = 2;
    auto mesh = createBeamMesh(v0, vEnd, nbBeams, 2);
    mesh->display(false);

            // siconos::mechanics::fem::writeMeshforPython(mesh);

    int bulk_material_tag = 1;
    int boundary_condition_tag = 2;
    int applied_force_tag = 3;

            // std::shared_ptr<Material> mat1 = std::make_shared<Material>(1, 8*36/5.,
            // 1/5.); // material for  triangle_felippa.msh
    double density = 7800.;
    auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(
        density, 210e9, 1 / 3);
    // auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(
    //     1.0, 1.0, 1 / 3);
    std::map<unsigned int, std::shared_ptr<siconos::mechanics::fem::Material>>
        materials = {{bulk_material_tag, mat1},{boundary_condition_tag,mat1},{applied_force_tag,mat1}};
    auto start = std::chrono::system_clock::now();
    std::cout << "Creating beamTIDS... " << std::endl;
    auto beam = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(
        mesh, materials, siconos::algebra::UblasType::SPARSE);
    // auto beam = std::make_shared<siconos::mechanics::fem::FiniteElementLinearTIDS>(
    //     mesh, materials, siconos::algebra::UblasType::SPARSE);
    std::cout << "beamTIDS created!" << std::endl;
    auto femodel = beam->FEModel();
    /*------------------------------------------------- Applied forces  */
                                                                            //    (*nodal_forces)(0) = 1e6;
    std::cout << "Applying forces :" << std::endl;

    auto nodal_forces = std::make_shared<Vector>(3);
    nodal_forces->zero();
    (*nodal_forces)(0) = -2e9;
    beam->applyNodalForces(applied_force_tag, nodal_forces);
    std::cout << "forces applied!" << std::endl;
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
    double T = 1e-02;    // final computation time
    double h = 1e-05;    // time step
    double theta = 1.0;  // theta for MoreauJeanOSI integrator
    int N = ceil((T - t0) / h);  // Number of time steps


    auto beamNSDS =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

            // add the dynamical system in the non smooth dynamical system
    beamNSDS->insertDynamicalSystem(beam);

    auto OSI = std::make_shared<siconos::integrators::MoreauJeanGOSI>(theta);
    // auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
    OSI->setIsWSymmetricDefinitePositive(true);

            // -- (2) Time discretisation --

    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::GlobalFrictionContact>(2, SICONOS_FRICTION_2D_NSGS);
    // auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();


            // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(beamNSDS, t,OSI,osnspb);
    // s->insertIntegrator(OSI);

    std::string SOFAfilename = "beam.state";
    int k = 1;
    start = std::chrono::system_clock::now();
    unsigned int outputSize = 4;
    Matrix dataPlot(N + 1, outputSize);
    auto q = beam->q();
    auto v = beam->velocity();

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
    dataPlot(0, 1) = (*q)(beam->dimension() - 3);
    dataPlot(0, 2) = (*q)(beam->dimension() - 2);
    dataPlot(0, 3) = (*q)(beam->dimension() - 1);

    std::cout << "About to start simulation" << std::endl;
    while (s->hasNextEvent()) {
      s->computeOneStep();
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(beam->dimension() - 3);
      dataPlot(k, 2) = (*q)(beam->dimension() - 2);
      dataPlot(k, 3) = (*q)(beam->dimension() - 1);
      vcnt = 0;
      // for (auto v : mesh->vertices()) {
      //   std::cout << v->num() << std::endl;

              //   pos(k,(vcnt)*2) = v->x() + (*q)((vcnt)*3);
              //   pos(k,(vcnt)*2+1) = v->y() + (*q)((vcnt)*3+1);
              //   vcnt++;
              // }
      writeBeamPositionforSOFA(mesh, femodel, q, SOFAfilename,k*0.01);

      s->nextStep();
      k++;

      siconos::tools::progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                       .count();


    siconos::algebra::io::write("bernoulliBeam-GOSI.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "bernoulliBeam.ref", eps)) >= eps)
      return 1;

    return 0;

  } catch (...) {
    std::cerr << "Exception caught in bernoulliBeam.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
