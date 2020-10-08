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

#include "./src/gmsh_io.hpp"

#include "FiniteElementLinearTIDS.hpp"
#include "Material.hpp"
//#include "FiniteElementModel.hpp"
//#include "FemFwd.hpp"

using namespace std;


static Mesh* createMesh2x1()
{
  MVertex * v1 = new MVertex(1, 0.,0.,0.);
  MVertex * v2 = new MVertex(2, 1.,0.,0.);
  MVertex * v3 = new MVertex(3, 0.,1.,0.);
  MVertex * v4 = new MVertex(4, 1.,1.,0.);

  std::vector<MVertex *> vertices = {v1, v2, v3, v4};


  std::vector<MVertex *> vertices1 = {v1, v2, v3};
  MElement * e1 = new MElement(1, 2, vertices1);

  std::vector<MVertex *> vertices2 = {v2, v4, v3};
  MElement * e2 = new MElement(2, 2, vertices2);

  std::vector<MElement *> elements = {e1, e2};

  return new Mesh(2, vertices, elements);
}
static Mesh* createMeshnxm(int n, int m, double Lx, double Ly)
{

  double lx = Lx/n;
  double ly = Ly/m;


  std::vector<MVertex *> vertices;

  vertices.resize((n+1)*(m+1));

  for (int i =0;  i < n+1; i++)
  {
    for (int j =0;  j < m+1; j++)
    {
      vertices[i+j*(n+1)] = new MVertex(i+j*(n+1), i*lx , j*ly, 0.);
    }
  }
  std::vector<MElement *> elements;
  //elements.resize(2*n*m);
  int element_cnt=0;
  for (int i =0;  i < n; i++)
  {
    for (int j =0;  j < m; j++)
    {
      std::vector<MVertex *> vertices_e_1 = {vertices[i+j*(n+1)], vertices[i+1+(j)*(n+1)], vertices[i+(j+1)*(n+1)]};
      elements.push_back(new MElement(element_cnt++, 2, vertices_e_1));

      std::vector<MVertex *> vertices_e_2 = {vertices[i+1+(j)*(n+1)], vertices[i+1+(j+1)*(n+1)], vertices[i+(j+1)*(n+1)]};
      elements.push_back(new MElement(element_cnt++, 2, vertices_e_2));
    }
  }
  return new Mesh(2, vertices, elements);
}
static Mesh* createMeshFromGMSH(string gmsh_filename)
{




  int *element_node;
  int element_num;
  int element_order;

  int m;
  int node_num;
  double *node_x;

  cout << "\n";
  cout << "  Read data from a file.\n";
//
//  Get the data size.
//
  gmsh_size_read ( gmsh_filename, node_num, m, element_num,
    element_order );
//
//  Print the sizes.
//
  cout << "\n";
  cout << "  Node data read from file \"" << gmsh_filename << "\"\n";
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";
  cout << "  Spatial dimension = " << m << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Element order = " << element_order << "\n";

//
//  Allocate memory.
//
  node_x = ( double * ) malloc ( m * node_num * sizeof ( double ) );
  element_node = ( int * )
    malloc ( element_order * element_num * sizeof ( int ) );
  std::vector<MVertex *> vertices;
  //vertices.resize(node_num);
  std::vector<MElement *> elements;
  //elements.resize(elements);
//
//  Get the data.
//
  gmsh_data_read ( gmsh_filename, m, node_num, node_x,
    element_order, element_num, element_node );

  for(int v = 0; v <node_num; v++)
  {
    vertices.push_back(new MVertex(v, node_x[0+v*m] , node_x[1+v*m], 0.));
  }
  unsigned int element_cnt=0;
  for(int e = 0; e <element_num; e++)
  {
    std::vector<MVertex *> vertices_e ;
    for ( int i = 0; i < element_order; i++ )
    {
      int num_v = element_node[i+e *element_order] - 1;
      vertices_e.push_back(vertices[num_v]);
    }

    elements.push_back(new MElement(element_cnt++, 2, vertices_e));
  }



//
//  Print some of the data.
//
  r8mat_transpose_print_some ( m, node_num, node_x,
    1, 1, m, 10, "  Coordinates for first 10 nodes:" );

  i4mat_transpose_print_some ( element_order, element_num, element_node,
    1, 1, element_order, 10, "  Connectivity for first 10 elements:" );
  return new Mesh(m, vertices, elements);

}
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
  fclose(foutput);
}

static void  outputDisplacementforPython(SP::Mesh  mesh, SP::FiniteElementModel femodel, SP::SiconosVector x)
{
  FILE * foutput = fopen("displacement.py", "a");
  fprintf(foutput, "x.append(np.array([");

  for (MVertex * v : mesh->vertices())
  {
    SP::FENode n = femodel->vertexToNode(v);
    unsigned int idx= (*n->dofIndex())[0];
    double value =(*x)(idx);
    fprintf(foutput, "%e,", value) ;

  }
  fprintf(foutput, "]))\n") ;
  fprintf(foutput, "\n");

  fprintf(foutput, "y.append(np.array([");
  for (MVertex * v : mesh->vertices())
  {
    SP::FENode n = femodel->vertexToNode(v);
    unsigned int idx= (*n->dofIndex())[1];
    double value =(*x)(idx);
    fprintf(foutput, "%e,", value) ;

  }
  fprintf(foutput, "]))\n") ;
  fprintf(foutput, "\n");
  fclose(foutput);
}


int main(int argc, char* argv[])
{
  
  double Ly= 0.3;
//  SP::Mesh mesh (createMesh2x1());
//  SP::Mesh mesh (createMeshnxm(50, 15 , 3., Ly));
  Ly =10.0;
  string gmsh_filename = "./step_2d.msh";
//   string gmsh_filename = "./square_v4.msh";
  SP::Mesh mesh (createMeshFromGMSH(gmsh_filename));
  mesh->display(true);
  outputMeshforPython(mesh);


  SP::Material material(new Material(1, 100000, 0.0));

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

    for(SP::FENode n : femodel->nodes())
    {
      //std::cout << "node number : " << n->num() << " " << n->x() << " " << n->y() <<  std::endl;
      if (fabs(n->y()-Ly) <= 1e-16 and fabs(n->x()) >= 1e-16)
      {
        //std::cout << "node number : " << n->num() << " " << n->y() <<  std::endl;
        unsigned int idx_y = (*n->dofIndex())[1];
        (*forces)(idx_y) = -10000.;
      }
    }
    //(*forces)(FEsolid->dimension()-1) = -100.;
    FEsolid->setFExtPtr(forces);


    /*------------------------------------------------- Boundary Conditions  */
    /* This part should be hidden in a new BC function for a node number
     * and a dof index. */
    SP::IndexInt bdIndex(new IndexInt(0));
    for(SP::FENode n : femodel->nodes())
    {
      if (fabs(n->x()) <= 1e-16)
      {
        //std::cout << "node number : " << n->num() << " " << n->x() <<  std::endl;
        unsigned int idx_x = (*n->dofIndex())[0];
        bdIndex->push_back(idx_x);
        unsigned int idx_y = (*n->dofIndex())[1];
        bdIndex->push_back(idx_y);
      }
    }
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
    double T = 1e-01;                  // final computation time
    double h = 1e-03;                // time step
    double theta = 1.0;              // theta for MoreauJeanOSI integrator

    SP::NonSmoothDynamicalSystem solid(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);

    /*------------------------------------------------- Contact Conditions  */
    double e =0.0;
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::SiconosVector  initial_gap(new SiconosVector(1, Ly*0.1));
    for(SP::FENode n : femodel->nodes())
    {
      if (fabs(n->y()) <= 1e-16 and fabs(n->x()) >= 1e-16)
      {
        //std::cout << "node number : " << n->num() << " " << n->y() <<  std::endl;
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
    ioMatrix::write("T3.dat", "ascii", dataPlot, "noDim");
    // double error=0.0, eps=1e-12;
    // if((error=ioMatrix::compareRefFile(dataPlot, "BouncingBallTS.ref", eps)) >= 0.0
    //     && error > eps)
    //   return 1;






  }
  catch(SiconosException& e)
  {
    cerr << e.report() << endl;
    return 1;

  }
  catch(...)
  {
    cerr << "Exception caught in T3.cpp" << endl;
    return 1;

  }




}
