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

//#include "cube_fine.cpp"
//#include "beam.cpp"

static void  outputMeshforPython(SP::Mesh  mesh)
{
  FILE * foutput = fopen("mesh.py", "w");
  fprintf(foutput, "coord=[]\n");
  for (MVertex * v : mesh->vertices())
  {
    fprintf(foutput, "coord.append([%e, %e])\n", v->x(), v->y());
  }
  fprintf(foutput, "triangle=[]\n");
  for (MElement * e : mesh->elements())
  {
    fprintf(foutput, "triangle.append([");
    for (MVertex * v : e->vertices())
    {
      fprintf(foutput, "%zu, ", v->num());
    }
    fprintf(foutput, "])\n");
  }
  fclose(foutput);
}

static void  preOutputDisplacementforPython()
{
  FILE * foutput = fopen("displacement.py", "w");
  fprintf(foutput, "import numpy as np\nx=[]\n");
  fprintf(foutput, "import numpy as np\ny=[]\n");
  fprintf(foutput, "import numpy as np\nz=[]\n");
  fclose(foutput);
}

static void  outputDisplacementforPython(SP::Mesh  mesh, SP::FiniteElementModel femodel, SP::SiconosVector x)
{
  FILE * foutput = fopen("displacement.py", "a");
  fprintf(foutput, "x.append(np.array([");

  for (MVertex * v : mesh->vertices())
  {
    SP::FENode n = femodel->vertexToNode(v);
    double value = 0.0;
    if (n)
    {
      unsigned int idx= (*n->dofIndex())[0];
      value =(*x)(idx);
    }
    fprintf(foutput, "%e,", value) ;

  }
  fprintf(foutput, "]))\n") ;
  
  fprintf(foutput, "\n");


  fprintf(foutput, "y.append(np.array([");
  for (MVertex * v : mesh->vertices())
  {
    SP::FENode n = femodel->vertexToNode(v);
    double value = 0.0;
    if (n)
    {
      unsigned int idx= (*n->dofIndex())[1];
      value =(*x)(idx);
    }
    fprintf(foutput, "%e,", value) ;

  }
  fprintf(foutput, "]))\n") ;
  fprintf(foutput, "\n");


  fprintf(foutput, "z.append(np.array([");
  for (MVertex * v : mesh->vertices())
  {
    SP::FENode n = femodel->vertexToNode(v);
    double value = 0.0;
    if (n)
    {
      unsigned int idx= (*n->dofIndex())[2];
      value =(*x)(idx);
    }
    fprintf(foutput, "%e,", value) ;

  }
  fprintf(foutput, "]))\n") ;
  fprintf(foutput, "\n");

  
  fclose(foutput);
}


int main(int argc, char* argv[])
{
  
  double Lz= 1.0;
  //SP::Mesh mesh (createMesh());
  SP::Mesh mesh (createMeshFromGMSH2("cube.msh2"));
//mesh->display(false);
  outputMeshforPython(mesh);

  SP::Material material(new Material(2500, 100000, 0.0));

  try{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    SP::FiniteElementLinearTIDS FEsolid (new FiniteElementLinearTIDS(mesh, material, Siconos::SPARSE));
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
      (end-start).count();
    cout << "Assembly time : " << elapsed << " ms" << endl;
    std::cout << " " << std::endl;
    //FEsolid->display(true);

    SP::FiniteElementModel femodel = FEsolid->FEModel();


    /*------------------------------------------------- Applied forces  */
    SP::SiconosVector forces(new SiconosVector(FEsolid->dimension()));
    forces->zero();
    std::cout << "node number applied forces: [";
    for(SP::FENode n : femodel->nodes())
    {
      //std::cout << "node number : " << n->num() << " " << n->x() << " " << n->y() <<  std::endl;
      if (fabs(n->z()-Lz) <= 1e-16 and fabs(n->x()) >= 1e-16)
      {
        std::cout   << " " << n->num();
        unsigned int idx_z = (*n->dofIndex())[2];
        (*forces)(idx_z) = -10;
      }
    }
    std::cout << " ]"  <<  std::endl;
    FEsolid->setFExtPtr(forces);


    /*------------------------------------------------- Boundary Conditions  */
    /* This part should be hidden in a new BC function for a node number
     * and a dof index. */
    SP::IndexInt bdIndex(new IndexInt(0));

    std::cout << "Boundary conditions node number : [ ";
    for(SP::FENode n : femodel->nodes())
    {
      if (fabs(n->x()) <= 1e-16)
      {
        std::cout  << n->num() << " " ;
        unsigned int idx_x = (*n->dofIndex())[0];
        bdIndex->push_back(idx_x);
        unsigned int idx_y = (*n->dofIndex())[1];
        bdIndex->push_back(idx_y);
        unsigned int idx_z = (*n->dofIndex())[2];
        bdIndex->push_back(idx_z);
      }
    }
    std::cout  <<  "] " <<  std::endl;
    SP::SiconosVector bdPrescribedVelocity(new SiconosVector(bdIndex->size()));
    for (int i=0; i < bdIndex->size() ; i++)
    {
      bdPrescribedVelocity->setValue(i,0.0);
    }
    SP::BoundaryCondition bd (new BoundaryCondition(bdIndex,bdPrescribedVelocity));
    FEsolid->setBoundaryConditions(bd);


    // -------------
    // --- Model ---
    // -------------
    double t0 = 0;                   // initial computation time
    double T = 100;                  // final computation time
    double h = 1e-00;                // time step
    double theta = 0.5;              // theta for MoreauJeanOSI integrator

    SP::NonSmoothDynamicalSystem solid(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);

    /*------------------------------------------------- Contact Conditions  */
    double e =0.0;
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::SiconosVector  initial_gap(new SiconosVector(1, Lz*0.05));
    std::cout << "contact node number : [ "  ;
    for(SP::FENode n : femodel->nodes())
    {
      if (fabs(n->z()) <= 1e-16 and fabs(n->x()) >= 8.)
      {
        std::cout << " " << n->num() ;
        unsigned int idx_z = (*n->dofIndex())[2];
        SP::SimpleMatrix H(new SimpleMatrix(1, FEsolid->dimension()));
        (*H)(0, idx_z) = 1.0;
        SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
        SP::Relation relation(new LagrangianLinearTIR(H, initial_gap));
        SP::Interaction inter(new Interaction(nslaw, relation));
        // link the interaction and the dynamical system
        solid->link(inter, FEsolid);
      }
    }
    std::cout << "]"<< std::endl;

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
    preOutputDisplacementforPython();



    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    start = std::chrono::system_clock::now();
    while(s->hasNextEvent())
    {
      s->computeOneStep();
      //osnspb->display();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      //std::cout << (*q)(0) << std::endl;
      dataPlot(k, 1) = (*q)(FEsolid->dimension()-1);
      dataPlot(k, 2) = (*v)(FEsolid->dimension()-1);
      dataPlot(k, 3) = (*p)(0);
      outputDisplacementforPython(mesh, femodel, q);
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
    ioMatrix::write("TH4.dat", "ascii", dataPlot, "noDim");
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
