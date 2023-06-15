
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
//-----------------------------------------------------------------------
//
//  CircuitRLCD  : sample of an electrical circuit involving :
//  - a linear dynamical system consisting of an LC oscillator (1 microF , 10 mH)
//  - a non smooth system (a 1000 Ohm resistor in series with a diode) in parallel
//    with the oscillator
//
//  Expected behavior :
//  The initial state of the oscillator provides an initial energy.
//  The period is 2 Pi sqrt(LC) ~ 0,628 ms.
//  A positive voltage across the capacitor allows current to flow
//  through the resistor-diode branch , resulting in an energy loss :
//  the oscillation damps.
//
//  State variables :
//  - the voltage across the capacitor (or inductor)
//  - the current through the inductor
//
//  Since there is only one dynamical system, the interaction is defined by :
//  - a complementarity law between diode current and voltage where y stands
//    for the reverse voltage across the diode and lambda stands for the
//    the diode current
//  - a linear time invariant relation between the state variables and
//    y and lambda (derived from Kirchhoff laws)
//
//-----------------------------------------------------------------------

#include "SiconosKernel.hpp"
#include <chrono>
#include <stdio.h>
#include <stdlib.h>

#include "MeshUtils.hpp"
#include "FiniteElementLinearTIDS.hpp"
#include "SiconosAlgebraProd.hpp"


using namespace std;
//using namespace siconos::mechanics::fem::native;


int main(int argc, char* argv[])
{
    double t0 = 0.0;
    double T = 500e-1;        // Total simulation time
    double h_step = 1.0e-2;// Time step

    double m = 1; // mass
    double k = 10; // stiffness of the spring
    double b = 1; // B matrix, trivial here

    double vInit = 0.3;    // initial velocity
    double xInit = 0.0;    // initial position
    double sigmaInit = k*xInit;    // initial stress
    double sigmaMax = 1.5;

    double Ly= 1.0;
    //  std::shared_ptr<Mesh> mesh = create2dMesh2x1();
    //  std::shared_ptr<Mesh> mesh = create2dMeshnxm(50, 15 , 3., Ly);
    //string gmsh_filename = "./mesh_data/triangle_felippa.msh";
//    string gmsh_filename = "./mesh_data/triangle_reference.msh";
//    string gmsh_filename = "./mesh_data/triangle.msh";
    string gmsh_filename = "./mesh_data/square_6.msh";
//      string gmsh_filename = "./mesh_data/square_200.msh";
    //string gmsh_filename = "./mesh_data/square_2720.msh";

    std::shared_ptr<siconos::mechanics::fem::native::Mesh> mesh(siconos::mechanics::fem::native::createMeshFromGMSH2(gmsh_filename));
    //mesh->display(false);

    writeMeshforPython(mesh);

    int bulk_material_tag = 1;
    int boundary_condition_tag = 2;
    int applied_force_tag = 3;

    //std::shared_ptr<Material> mat1 = std::make_shared<Material>(1, 8*36/5., 1/5.); // material for  triangle_felippa.msh
    double density = 7000;
    std::shared_ptr<Material> mat1 = std::make_shared<Material>(density, 1.0e3, 0.3);
    std::map<unsigned int, std::shared_ptr<Material> > materials = {{bulk_material_tag, mat1}};

    std::shared_ptr<siconos::mechanics::fem::native::FiniteElementLinearTIDS> FEsolid  = std::make_shared<siconos::mechanics::fem::native::FiniteElementLinearTIDS>(mesh, materials, Siconos::SPARSE);
    std::shared_ptr<siconos::mechanics::fem::native::FiniteElementModel> femodel = FEsolid->FEModel();


    string Modeltitle = "SpringMass";

    try
    {
        int nElements = (femodel->elements()).size();
        double ndofs = (mesh->vertices()).size()*mesh->dim();
        double dim = 2*ndofs + 3*nElements;
        // --- Dynamical system specification ---
        SP::SiconosVector init_state(new SiconosVector(dim));
        SP::SiconosVector Fext(new SiconosVector(dim));
        Fext->setValue(3, 100.0);
        Fext->setValue(9, 100.0);
        SP::SimpleMatrix mass(new SimpleMatrix(ndofs, ndofs));
        SP::SimpleMatrix stiffness(new SimpleMatrix(ndofs, ndofs));
        SP::SimpleMatrix Bfem(new SimpleMatrix(3, ndofs));
        femodel->computeMassMatrix(mass, materials);
        std::shared_ptr<SimpleMatrix> D = std::make_shared<SimpleMatrix>(3,3);
        double E = mat1->elasticYoungModulus();
        double nu =  mat1->poissonCoefficient();

        double coef = E/((1+nu)*(1-2.*nu));
        (*D)(0,0) = coef*(1.-nu);
        (*D)(0,1) = coef*nu;
        (*D)(0,2) = 0.0;

        (*D)(1,0) = (*D)(0,1);
        (*D)(1,1) = (*D)(0,0);
        (*D)(1,2) = 0.0;

        (*D)(2,0) = 0.0;
        (*D)(2,1) = 0.0;
        (*D)(2,2) = 0.5*coef*(1.0 - 2* nu);

        //    double coef = E/(1-nu*nu);
        //    (*D)(0,0) = coef;
        //    (*D)(0,1) = coef*nu;
        //    (*D)(0,2) = 0.0;

        //    (*D)(1,0) = (*D)(0,1);
        //    (*D)(1,1) = (*D)(0,0);
        //    (*D)(1,2) = 0.0;

        //    (*D)(2,0) = 0.0;
        //    (*D)(2,1) = 0.0;
        //    (*D)(2,2) = 0.5*coef*(1.0 - nu);

        femodel->computeStiffnessMatrix(stiffness, materials);

        std::cout << "nElements:" << nElements << std::endl;
        SP::SimpleMatrix BigB(new SimpleMatrix(3*nElements, ndofs));
        std::shared_ptr<SimpleMatrix> DBigB = std::make_shared<SimpleMatrix>(3*nElements,ndofs);
        std::shared_ptr<SimpleMatrix> DBfem = std::make_shared<SimpleMatrix>(3,ndofs);
        int elem_cnt = 0;
        for(std::shared_ptr<siconos::mechanics::fem::native::FElement> fe : femodel->elements())
        {
                    femodel->computeB_Matrix_direct(*fe, *Bfem);
                    prod(*D, *Bfem, *DBfem, true);
                    femodel->AssembleElementary_B_Matrix(BigB,*Bfem,*fe, elem_cnt);
                    femodel->AssembleElementary_B_Matrix(DBigB,*DBfem,*fe, elem_cnt);
                    elem_cnt++;
        }
        SP::SimpleMatrix M(new SimpleMatrix(dim, dim));
        for (int i=0;i<ndofs;i++){
            for (int j=0;j<ndofs;j++){
                M->setValue(i, j, mass->getValue(i,j));
            }
        }
        for (int i=0;i<ndofs;i++){
            M->setValue(i+ndofs, i+ndofs, 1.0);
        }

        for (int i=0;i<3*nElements;i++){
            M->setValue(i+2*ndofs, i+2*ndofs, 1.0);
        }

        std::shared_ptr<SimpleMatrix> BT = std::make_shared<SimpleMatrix>(ndofs,3*nElements);
        BT->trans(*BigB);
        std::shared_ptr<SimpleMatrix> BTDB = std::make_shared<SimpleMatrix>(ndofs,ndofs);

        prod(*BT, *DBigB, *BTDB, true);

        SP::SiconosVector velocity(new SiconosVector(ndofs));
        SP::SiconosVector DBv(new SiconosVector(3));
        SP::SiconosVector DBu(new SiconosVector(3));

        SP::SimpleMatrix A(new SimpleMatrix(dim, dim));
        for (int i=0;i<ndofs;i++){
            for (int j=0;j<3*nElements;j++){
                A->setValue(i, j+2*ndofs, -BigB->getValue(j,i));
            }
        }
        for (int i=0;i<3*nElements;i++){
            for (int j=0;j<ndofs;j++){
                A->setValue(i+2*ndofs, j, DBigB->getValue(i,j));
            }
        }
        for (int i=0;i<ndofs;i++){
            A->setValue(i+ndofs, i, 1.0);
        }

        SP::FirstOrderLinearTIDS triangle(new FirstOrderLinearTIDS(init_state, A, Fext));
        triangle->setMPtr(M);
        triangle->display();
        std::cout << "q and velocity from FEmodel:" << std::endl;
        std::shared_ptr<SiconosVector> q = FEsolid->q();
        std::shared_ptr<SiconosVector> velo = FEsolid->velocity();
        for (int i=0;i<6;i++)
            std::cout << (*q)(i) << std::endl;
        std::cout << "Velocity from FEmodel:" << std::endl;
        for (int i=0;i<6;i++)
            std::cout << (*velo)(i) << std::endl;
        // --- Interaction between linear system and non smooth system ---
        int nbConstraints = 10;
        SP::SimpleMatrix C(new SimpleMatrix(nbConstraints, dim));
        // Position Constraints
        C->setValue(0, ndofs, -1.0);
        C->setValue(1, ndofs, 1.0);
        C->setValue(2, ndofs+1, -1.0);
        C->setValue(3, ndofs+1, 1.0);
        C->setValue(4, ndofs+6, -1.0);
        C->setValue(5, ndofs+6, 1.0);
        C->setValue(6, ndofs+7, -1.0);
        C->setValue(7, ndofs+7, 1.0);
        // Stress Constraints
        C->setValue(8, 2*ndofs, -1.0);
        C->setValue(9, 2*ndofs, 1.0);

        SP::SimpleMatrix B(new SimpleMatrix(dim, nbConstraints));
        B->setValue(0, 0, -1.0);
        B->setValue(0, 1, 1.0);
        B->setValue(0+1, 2, -1.0);
        B->setValue(0+1, 3, 1.0);
        B->setValue(0+6, 4, -1.0);
        B->setValue(0+6, 5, 1.0);
        B->setValue(0+7, 6, -1.0);
        B->setValue(0+7, 7, 1.0);
        // Stress Constraints
        B->setValue(2*ndofs, 8, -1.0);
        B->setValue(2*ndofs, 9, 1.0);
        SP::SiconosVector e(new SiconosVector(nbConstraints));
        e->setValue(0, 0.00001);
        e->setValue(1, 0.00001);
        e->setValue(2, 0.00001);
        e->setValue(3, 0.00001);
        e->setValue(4, 0.00001);
        e->setValue(5, 0.00001);
        e->setValue(6, 0.00001);
        e->setValue(7, 0.00001);
        // Stress Constraints
        e->setValue(8, 100);
        e->setValue(9, 100);
        SP::FirstOrderLinearTIR LTIRspring(new FirstOrderLinearTIR(C, B));
        LTIRspring->setePtr(e);
        SP::NonSmoothLaw NSLaw(new ComplementarityConditionNSL(nbConstraints));

        SP::Interaction InterTriangle(new Interaction(NSLaw, LTIRspring));
        InterTriangle->display();
        // --- Model creation ---
        SP::NonSmoothDynamicalSystem springDS(new NonSmoothDynamicalSystem(t0, T));
        springDS->setTitle(Modeltitle);
        // add the dynamical system in the non smooth dynamical system
        springDS->insertDynamicalSystem(triangle);

        // link the interaction and the dynamical system
        springDS->link(InterTriangle, triangle);

        InterTriangle->computeOutput(t0,0);
        InterTriangle->computeInput(t0,0);


        springDS->display();

        // ------------------
        // --- Simulation ---
        // ------------------
        double theta = 0.5000000000001;
        //    double theta = 1.0;

        // -- (1) OneStepIntegrators --
        SP::EulerMoreauOSI OSI(new EulerMoreauOSI(theta));

        // -- (2) Time discretisation --
        SP::TimeDiscretisation TiDis(new TimeDiscretisation(t0, h_step));
        // --- (3) one step non smooth problem
        SP::LCP LCP_triangle(new LCP());
        //    SP::Relay LCP_spring(new Relay());

        // -- (4) Simulation setup with (1) (2) (3)
        SP::TimeStepping springMassTS(new TimeStepping(springDS, TiDis,OSI ,LCP_triangle));
        //    SP::TimeStepping springMassTS(new TimeStepping(springDS, TiDis));
        //    springMassTS->insertIntegrator(OSI);
        double h = springMassTS->timeStep();
        int N = ceil((T - t0) / h); // Number of time steps
        int k = 0;

        // --- Get the values to be plotted ---
        // -> saved in a matrix dataPlot
        SimpleMatrix dataPlot(N, 22);

        // For the initial time step:
        std::cout << "q and velocity initially:" << std::endl;
        for (int i=6;i<12;i++)
            std::cout << (*triangle->x())(i) << std::endl;
        std::cout << "Velocity from FEmodel:" << std::endl;
        for (int i=0;i<6;i++)
            std::cout << (*triangle->x())(i) << std::endl;

        // time
        dataPlot(k, 0) = springMassTS->nextTime();
        dataPlot(k, 1) = (*triangle->x())(2*ndofs);
        dataPlot(k, 2) = (*triangle->x())(2*ndofs+1);
        dataPlot(k, 3) = (*triangle->x())(2*ndofs+2);
        dataPlot(k, 4) = (*triangle->x())(2);
        dataPlot(k, 5) = (*triangle->x())(3);
        dataPlot(k, 6) = (*triangle->x())(ndofs+2);
        dataPlot(k, 7) = (*triangle->x())(ndofs+3);
        dataPlot(k, 8) = (InterTriangle->getLambda(0))(8);
        dataPlot(k, 9) = (InterTriangle->getLambda(0))(9);

        double x,v,sigma,plasticRate,epElastic,ekinetic,plasticDeformation=0.0,plasticDissipation=0.0;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        // --- Time loop  ---
        SP::SiconosVector displacement(new SiconosVector(ndofs));
        std::string filename = "triangle.state";
        FILE * foutput = fopen(filename.c_str(), "w");
        fclose(foutput);
        for(k = 1 ; k < N ; ++k)
        {
            // solve ...
            springMassTS->computeOneStep();
            std::shared_ptr<SiconosVector> q = FEsolid->q();
            std::shared_ptr<SiconosVector> velo = FEsolid->velocity();
            for (int i=0;i<ndofs;i++)
                displacement->setValue(i,(*triangle->x())(i+ndofs));
            siconos::mechanics::fem::native::writePositionforSOFA(mesh,femodel, displacement, filename,  springMassTS->nextTime());
            // --- Get values to be plotted ---
            // time
            dataPlot(k, 0) = springMassTS->nextTime();
            dataPlot(k, 1) = (*triangle->x())(2*ndofs);
            dataPlot(k, 2) = (*triangle->x())(2*ndofs+1);
            dataPlot(k, 3) = (*triangle->x())(2*ndofs+2);
            dataPlot(k, 4) = (*triangle->x())(2);
            dataPlot(k, 5) = (*triangle->x())(3);
            dataPlot(k, 6) = (*triangle->x())(ndofs+2);
            dataPlot(k, 7) = (*triangle->x())(ndofs+3);
            dataPlot(k, 8) = (InterTriangle->getLambda(0))(8);
            dataPlot(k, 9) = (InterTriangle->getLambda(0))(9);

            // transfer of state i+1 into state i and time incrementation
            springMassTS->nextStep();

        }
        // Number of time iterations
        cout << "Number of iterations done: " << k - 1 << endl;
        cout << "Computation Time " << endl;
        end = std::chrono::system_clock::now();
        int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                (end-start).count();
        cout << "Computation time : " << elapsed << " ms" << endl;


        // dataPlot (ascii) output
        ioMatrix::write("plasticTriangle.dat", "ascii", dataPlot, "noDim");

        //    double error=0.0, eps=1e-12;
        //    if((error=ioMatrix::compareRefFile(dataPlot, "CircuitRLCD.ref", eps)) >= 0.0
        //        && error > eps)
        //      return 1;

    }


    // --- Exceptions handling ---
    catch(...)
    {
        Siconos::exception::process();
        return 1;
    }
}
