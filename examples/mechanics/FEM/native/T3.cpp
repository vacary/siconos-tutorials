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

#include "MeshUtils.hpp"

#include "FiniteElementLinearTIDS.hpp"
#include "Material.hpp"
//#include "FiniteElementModel.hpp"
//#include "FemFwd.hpp"

using namespace std;

int main(int argc, char* argv[])
{

  double Ly= 1.0;
//  SP::Mesh mesh (create2dMesh2x1());
//  SP::Mesh mesh (create2dMeshnxm(50, 15 , 3., Ly));
  Ly =1.0;
  //string gmsh_filename = "./mesh_data/triangle_felippa.msh";
  //string gmsh_filename = "./mesh_data/triangle_reference.msh";
  //string gmsh_filename = "./mesh_data/square_6.msh";
  string gmsh_filename = "./mesh_data/square_2720.msh";

  SP::Mesh mesh (createMeshFromGMSH2(gmsh_filename));
  mesh->display(false);

  writeMeshforPython(mesh);

  int bulk_material_tag = 1;
  int boundary_condition_tag = 2;
  int applied_force_tag = 3;

  //SP::Material mat1(new Material(1, 8*36/5., 1/5.)); // material for  triangle_felippa.msh
  double density = 7800.;
  SP::Material mat1(new Material(density, 210e9, 1/3.));
  std::map<unsigned int, SP::Material> materials = {{bulk_material_tag, mat1}};


  try{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    SP::FiniteElementLinearTIDS FEsolid (new FiniteElementLinearTIDS(mesh, materials, Siconos::SPARSE));
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
      (end-start).count();
    cout << "Assembly time : " << elapsed << " ms" << endl;
    std::cout << " " << std::endl;
    //FEsolid->display(true);

    SP::FiniteElementModel femodel = FEsolid->FEModel();
    // FEsolid->K()->display();
    // getchar();


    /*------------------------------------------------- Applied forces  */

    SP::SiconosVector nodal_forces(new SiconosVector(2));
    nodal_forces->zero();
    //(*nodal_forces)(0) = 1e6;
    (*nodal_forces)(1) = -1e7;
    FEsolid->applyNodalForces(applied_force_tag , nodal_forces);


    /*------------------------------------------------- Boundary Conditions  */
    /* This part should be hidden in a new BC function for a node number
     * and a dof index. */

    SP::IndexInt node_dof_index(new IndexInt(0));
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

    SP::NonSmoothDynamicalSystem solid(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);

    /*------------------------------------------------- Contact Conditions  */
    double e =0.0;
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::SiconosVector  initial_gap(new SiconosVector(1, Ly*5e-4));
    for(SP::FENode n : femodel->nodes())
    {
      if (fabs(n->y()) <= 1e-16 and fabs(n->x()) >= 1e-16)
      {
        std::cout << "contact node number : " << n->num() << " " << n->y() <<  std::endl;
        unsigned int idx_y = (*n->dofIndex())[1];
        SP::SimpleMatrix H(new SimpleMatrix(1, FEsolid->dimension()));
        (*H)(0, idx_y) = 1.0;
        SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
        SP::Relation relation(new LagrangianLinearTIR(H, initial_gap));
        SP::Interaction inter(new Interaction(nslaw, relation));
        // link the interaction and the dynamical system
        solid->link(inter, FEsolid);
      }
    }

    // // link the interaction and the dynamical system
    // bouncingBall->link(inter, FEsolid);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));
    OSI->setIsWSymmetricDefinitePositive(true);


    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(solid, t, OSI, osnspb));

    // =========================== End of model definition ===========================

    // ================================= Computation =================================


    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 5;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = FEsolid->q();
    SP::SiconosVector v = FEsolid->velocity();
    SP::SiconosVector p = FEsolid->p(1);
    //SP::SiconosVector lambda = inter->lambda(1);

    dataPlot(0, 0) = solid->t0();
    dataPlot(0, 1) = (*q)(FEsolid->dimension()-1);
    dataPlot(0, 2) = (*v)(FEsolid->dimension()-1);
    dataPlot(0, 3) = (*p)(0);
    //dataPlot(0, 4) = (*lambda)(0);

    std::string filename = prepareWriteDisplacementforPython("T3");
    writeDisplacementforPython(mesh, femodel, q, filename);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    start = std::chrono::system_clock::now();
    while(s->hasNextEvent())
    {
      s->computeOneStep();
      osnspb->display();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      //std::cout << (*q)(0) << std::endl;
      dataPlot(k, 1) = (*q)(FEsolid->dimension()-1);
      dataPlot(k, 2) = (*v)(FEsolid->dimension()-1);
      dataPlot(k, 3) = (*p)(0);

      if (k%1 == 0)
        writeDisplacementforPython(mesh, femodel, q, filename);
      //dataPlot(k, 4) = (*lambda)(0);
      s->nextStep();
      k++;
      progressBar((double)k/N);

    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
      (end-start).count();
    cout << endl <<  "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms" << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("T3.dat", "ascii", dataPlot, "noDim");
    // double error=0.0, eps=1e-12;
    // if((error=ioMatrix::compareRefFile(dataPlot, "BouncingBallTS.ref", eps)) >= 0.0
    //     && error > eps)
    //   return 1;






  }
  catch(...)
  {
   cerr << "Exception caught in TH4.cpp" << endl;
   Siconos::exception::process();
   return 1;
  }

}
