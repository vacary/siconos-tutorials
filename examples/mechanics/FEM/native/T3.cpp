/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*!\file T3.cpp
 */

#include <SiconosKernel.hpp>
#include <chrono>
#include <stdio.h>

#include "FENode.hpp"
#include "FiniteElementModel.hpp"
#include "MeshUtils.hpp"
#include "FiniteElementLinearTIDS.hpp"
#include "Material.hpp"

using namespace std;
using namespace siconos::mechanics::fem;

int main(int argc, char* argv[])
{

  double Ly= 1.0;
//  std::shared_ptr<Mesh> mesh = create2dMesh2x1();
//  std::shared_ptr<Mesh> mesh = create2dMeshnxm(50, 15 , 3., Ly);
  Ly =1.0;
  //string gmsh_filename = "./mesh_data/triangle_felippa.msh";
//  string gmsh_filename = "./mesh_data/triangle_reference.msh";
  string gmsh_filename = "./mesh_data/square_6.msh";
//  string gmsh_filename = "./mesh_data/square_200.msh";
  //string gmsh_filename = "./mesh_data/square_2720.msh";

  std::shared_ptr<Mesh> mesh(createMeshFromGMSH2(gmsh_filename));
  //mesh->display(false);

  writeMeshforPython(mesh);

  int bulk_material_tag = 1;
  int boundary_condition_tag = 2;
  int applied_force_tag = 3;

  //std::shared_ptr<Material> mat1 = std::make_shared<Material>(1, 8*36/5., 1/5.); // material for  triangle_felippa.msh
  double density = 7800.;
  std::shared_ptr<Material> mat1 = std::make_shared<Material>(density, 210e9, 1/3.);
  std::map<unsigned int, std::shared_ptr<Material> > materials = {{bulk_material_tag, mat1}};


  try
  {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    std::shared_ptr<FiniteElementLinearTIDS> FEsolid  = std::make_shared<FiniteElementLinearTIDS>(mesh, materials, siconos::algebra::UblasType::SPARSE);
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    cout << "Assembly time : " << elapsed << " ms" << endl;
    std::cout << " " << std::endl;
    //FEsolid->display(true);

    std::shared_ptr<FiniteElementModel> femodel = FEsolid->FEModel();
    // FEsolid->K()->display();
    // getchar();


    /*------------------------------------------------- Applied forces  */

    std::shared_ptr<siconos::algebra::SiconosVector> nodal_forces = std::make_shared<siconos::algebra::SiconosVector>(2);
    nodal_forces->zero();
//    (*nodal_forces)(0) = 1e6;
    (*nodal_forces)(1) = -1e7;
    FEsolid->applyNodalForces(applied_force_tag, nodal_forces);


    /*------------------------------------------------- Boundary Conditions  */
    /* This part should be hidden in a new BC function for a node number
     * and a dof index. */

    std::shared_ptr<std::vector<int>> node_dof_index = std::make_shared<std::vector<int>>(0);
    node_dof_index->push_back(0);
    node_dof_index->push_back(1);

    FEsolid->applyDirichletBoundaryConditions(boundary_condition_tag, node_dof_index);
    FEsolid->boundaryConditions()->display();

    // -------------
    // --- Model ---
    // -------------
    double t0 = 0;                   // initial computation time
    double T = 1e-02;                  // final computation time
    double h = 1e-05;                // time step
    double theta = 1.0;              // theta for MoreauJeanOSI integrator

    std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> solid = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    cout << "plop .. " << endl;
    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);

    /*------------------------------------------------- Contact Conditions  */
    double e =0.0;
    std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    std::shared_ptr<siconos::algebra::SiconosVector>  initial_gap = std::make_shared<siconos::algebra::SiconosVector>(1, Ly*5e-4);
    for(std::shared_ptr<FENode> n : femodel->nodes())
    {
      if(fabs(n->y()) <= 1e-16 and fabs(n->x()) >= 1e-16)
      {
        std::cout << "contact node number : " << n->num() << " " << n->y() <<  std::endl;
        unsigned int idx_y = (*n->dofIndex())[1];
        std::shared_ptr<siconos::algebra::SimpleMatrix> H = std::make_shared<siconos::algebra::SimpleMatrix>(1, FEsolid->dimension());
        (*H)(0, idx_y) = 1.0;
        std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
        std::shared_ptr<siconos::modeling::Relation> relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H, initial_gap);
        std::shared_ptr<siconos::modeling::Interaction> inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);
        solid->display();
        // link the interaction and the dynamical system
        solid->link(inter, FEsolid);
      }
    }
    cout << "plopdsaf .. " << endl;

    // // link the interaction and the dynamical system
    // bouncingBall->link(inter, FEsolid);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    std::shared_ptr<siconos::integrators::MoreauJeanOSI> OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
    OSI->setIsWSymmetricDefinitePositive(true);


    // -- (2) Time discretisation --
    std::shared_ptr<siconos::simulation::TimeDiscretisation> t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    std::shared_ptr<siconos::simulation::TimeStepping> s = std::make_shared<siconos::simulation::TimeStepping>(solid, t, OSI, osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    cout << "plosfdasfasfp .. " << endl;

    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 5;
    siconos::algebra::SimpleMatrix dataPlot(N + 1, outputSize);

    std::shared_ptr<siconos::algebra::SiconosVector> q = FEsolid->q();
    std::shared_ptr<siconos::algebra::SiconosVector> v = FEsolid->velocity();
    std::shared_ptr<siconos::algebra::SiconosVector> p = FEsolid->p(1);
//    std::shared_ptr<siconos::algebra::SiconosVector> lambda = inter->lambda(1);
    dataPlot(0, 0) = solid->t0();
    dataPlot(0, 1) = (*q)(FEsolid->dimension()-1);
    dataPlot(0, 2) = (*v)(FEsolid->dimension()-1);
//    dataPlot(0, 4) = (*lambda)(0);

    std::string filename = prepareWriteDisplacementforPython("T3");
    cout << "plop18 .. " << endl;
    writeDisplacementforPython(mesh, femodel, q, filename);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    start = std::chrono::system_clock::now();
    while(s->hasNextEvent())
    {
        cout << "k is ... : " << k << endl;

      s->computeOneStep();
      cout << "Out compute onestep ... : "  << endl;
      //osnspb->display();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      //std::cout << (*q)(0) << std::endl;
      dataPlot(k, 1) = (*q)(FEsolid->dimension()-1);
      dataPlot(k, 2) = (*v)(FEsolid->dimension()-1);
      dataPlot(k, 3) = (*p)(0);

      if(k%1 == 0)
        writeDisplacementforPython(mesh, femodel, q, filename);
      //dataPlot(k, 4) = (*lambda)(0);
      s->nextStep();
      k++;
      siconos::tools::progressBar((double)k/N);

    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
              (end-start).count();
    cout << endl <<  "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms" << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("T3.dat", dataPlot, siconos::algebra::io::ASCII_OUT, siconos::algebra::io::WriteType::nodim);
    double error=0.0, eps=1e-12;
    if((error=siconos::algebra::io::compareRefFile(dataPlot, "T3_square_200.ref", eps)) >= 0.0
        && error > eps)
      return 1;






  }
  catch(...)
  {
    cerr << "Exception caught in T3.cpp" << endl;
    siconos::exception::process();
    return 1;
  }

}
