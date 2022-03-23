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

#include "SiconosKernel.hpp"
#include <chrono>
#include "SiconosAlgebraProd.hpp"
//#define TS_PROJ
//#define TS_COMBINED
using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
#include "UserDefinedParameter.hpp"
//#include "UserDefinedParameter-ref.hpp"
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." <<endl<<endl;

    if(argc==3)
    {
      //printf("argv[0] %s\n", argv[0]);
      printf("h is set to %e\n", atof(argv[1]));
      printf("nElement is set to %i\n", atoi(argv[2]));
      h  = atof(argv[1]);
      nElement = atoi(argv[2]);
    }
    unsigned int nDof = (nElement +1);
    cout << "number of element = " << nElement <<  endl;
    cout << "number of dof = " << nDof <<  endl;

    double l = L/nElement; // length of an element

    cout << "bar length " << l <<  endl;

    SP::SiconosMatrix SparseMass(new SimpleMatrix(nDof,nDof,Siconos::SPARSE,nDof));
    SP::SiconosMatrix SparseStiffness(new SimpleMatrix(nDof,nDof,Siconos::SPARSE,3*nDof));


    SparseMass->setValue(0,0,1.0/3.0);
    SparseMass->setValue(0,1,1.0/6.0);
    SparseStiffness->setValue(0, 0, 1.0);
    SparseStiffness->setValue(0, 1, -1.0);

    for(unsigned int i = 1; i < nDof-1; i++)
    {
      SparseMass->setValue(i,i,2.0/3.0);
      SparseMass->setValue(i,i-1,1.0/6.0);
      SparseMass->setValue(i,i+1,1.0/6.0);


      SparseStiffness->setValue(i,i,2.0);
      SparseStiffness->setValue(i,i-1,-1.0);
      SparseStiffness->setValue(i,i+1,-1.0);
    }

    SparseMass->setValue(nDof-1,nDof-1,1.0/3.0);
    SparseMass->setValue(nDof-1,nDof-2,1.0/6.0);


    SparseStiffness->setValue(nDof-1, nDof-2,-1.0);
    SparseStiffness->setValue(nDof-1,nDof-1,1.0);

    *SparseMass  *= rho*S*l;
    *SparseStiffness  *= E*S/l;



//      SparseMass->display();
//      SparseStiffness->display();


    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof,position_init));
    SP::SiconosVector v0(new SiconosVector(nDof,velocity_init));

    // -- The dynamical system --
    SP::LagrangianLinearTIDS bar(new LagrangianLinearTIDS(q0,v0,SparseMass));

    // -- Set stiffness matrix (weight) --
    bar->setKPtr(SparseStiffness);


    // -- Set external forces (weight) --
    //SP::SiconosVector weight(new SiconosVector(nDof,-g*rho*S/l));
    SP::SiconosVector weight(new SiconosVector(nDof,0.0));
    bar->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.0;

    // Interaction bar-floor
    //
    SP::SimpleMatrix H(new SimpleMatrix(1,nDof));
    (*H)(0,0) = 1.0;

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianLinearTIR(H));

    SP::Interaction inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem impactingBar(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    impactingBar->insertDynamicalSystem(bar);

    // link the interaction and the dynamical system
    impactingBar->link(inter,bar);


    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
#ifdef TS_PROJ
    SP::MoreauJeanDirectProjectionOSI OSI(new MoreauJeanDirectProjectionOSI(theta));
    OSI->setDeactivateYPosThreshold(1e-05);
    OSI->setDeactivateYVelThreshold(0.0);
    OSI->setActivateYPosThreshold(1e-09);
    OSI->setActivateYVelThreshold(100.0);
#else
#ifdef TS_COMBINED
    SP::OneStepIntegrator OSI(new MoreauJeanCombinedProjectionOSI(theta));
#else
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta,0.0));
    //OSI->setConstraintActivationThreshold(1e-05);
#endif
#endif
    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0,h));

    // -- (3) one step non smooth problem
    SP::LCP osnspb(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
#ifdef TS_PROJ
    SP::MLCPProjectOnConstraints position(new MLCPProjectOnConstraints());
    SP::TimeSteppingDirectProjection s(new TimeSteppingDirectProjection(impactingBar, t,OSI, osnspb, position,0));
    s->setProjectionMaxIteration(10);
    s->setConstraintTolUnilateral(1e-10);
    s->setConstraintTol(1e-10);


#else
#ifdef TS_COMBINED
    SP::OneStepNSProblem position(new MLCPProjectOnConstraints(SICONOS_MLCP_ENUM));
    SP::TimeSteppingCombinedProjection s(new TimeSteppingCombinedProjection(impactingBar, t,OSI, osnspb, position,2));
    s->setProjectionMaxIteration(500);
    s->setConstraintTolUnilateral(1e-10);
    s->setConstraintTol(1e-10);
#else
    SP::TimeStepping s(new TimeStepping(impactingBar, t,OSI,osnspb));
#endif
#endif

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---


    int N = floor((T-t0)/h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 12;
    SimpleMatrix dataPlot(2,outputSize);

    SP::SiconosVector q = bar->q();
    SP::SiconosVector v = bar->velocity();
    SP::SiconosVector p = bar->p(1);
    SP::SiconosVector lambda = inter->lambda(1);

    dataPlot(0, 0) = impactingBar->t0();
    dataPlot(0,1) = (*q)(0);
    dataPlot(0,2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0,7) = (*q)(nDof-1);
    dataPlot(0,8) = (*v)(nDof-1);
    dataPlot(0,9) = (*q)((nDof)/2);
    dataPlot(0,10) = (*v)((nDof)/2);


    SP::SiconosVector tmp(new SiconosVector(nDof));

    prod(*SparseStiffness, *q, *tmp, true);
    double potentialEnergy = inner_prod(*q,   *tmp);
    prod(*SparseMass, *v, *tmp, true);
    double kineticEnergy = inner_prod(*v,*tmp);
    double impactEnergy = 0.0;

    dataPlot(0, 5) = potentialEnergy;
    dataPlot(0, 6) = kineticEnergy;
    dataPlot(0, 11) = impactEnergy;

//    std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
//    std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;



    // --- Time loop ---
    cout << "====> Start computation ... " <<endl<<endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;


    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    double max_impulse = 0.0;
    double W_00 =0.0;
    std::cout << "N :"  << N << std::endl;
//    while (s->nextTime() < T)
    while(k < 10)
    {
      // if (k%100 == 0)
      std::cout << "k :"  << k << std::endl;
      s->computeOneStep();

//       std::cout << "position"  << std::endl;
//       q->display();
//       std::cout << "velocity"  << std::endl;
//       v->display();

// --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k,1) = (*q)(0);
      dataPlot(k,2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0)/h;
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k,7) = (*q)(nDof-1);
      dataPlot(k,8) = (*v)(nDof-1);
      dataPlot(k,9) = (*q)((nDof)/2);
      dataPlot(k,10) = (*v)((nDof)/2);



      prod(*SparseStiffness, *q, *tmp, true);
      potentialEnergy = inner_prod(*q,   *tmp);
      prod(*SparseMass, *v, *tmp, true);
      kineticEnergy = inner_prod(*v,*tmp);

      dataPlot(k, 5) = potentialEnergy;
      dataPlot(k, 6) = kineticEnergy;

      double v_p_theta= theta*(*v)(0) + (1-theta) *  dataPlot(k-1,2);
      impactEnergy = v_p_theta * (*lambda)(0);
      dataPlot(k, 11) = impactEnergy;

//      std::cout << "q" << std::endl;
//       q->display();


//       std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
//       std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;
      if ((*lambda)(0) >0 )
      {
        max_impulse = fmax(max_impulse, (*lambda)(0));
        osnspb->display();
        SP::NumericsMatrix nm = osnspb->M()->numericsMatrix();
        W_00= NM_get_value(&(*nm),0,0);
        std::cout<< "W(0,0)" << NM_get_value(&(*nm),0,0) << std::endl;
        break;
      }


      s->nextStep();

      k++;
    }

    double totalImpactEnergy=0.0;
    for (int p =0; p < k; p++)
    {
      totalImpactEnergy = totalImpactEnergy + dataPlot(p, 11);
    }
    printf("totalImpactEnergy = %e\n", totalImpactEnergy);
    printf("max_impulse = %e\n", max_impulse);

    ofstream myfile;
    myfile.open ("percussion.txt", ios::out | ios::app );
    myfile  << nElement << " "
           << h << " "
           << totalImpactEnergy << " "
           << std::setprecision(24) << std::scientific 
           << max_impulse << " "
           << W_00 <<  endl;;
    myfile.close();
    exit(0);

    cout<<endl << "End of computation - Number of iterations done: "<<k-1<<endl;
    cout << "Computation Time " << endl;;
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    cout << "Computation time : " << elapsed << " ms" << endl;
    // --- Output files ---
    cout<<"====> Output file writing ..."<<endl;
    ioMatrix::write("ImpactingBar.dat", "ascii", dataPlot,"noDim");
    // cout << " Comparison with a reference file" << endl;
    // SimpleMatrix dataPlotRef(dataPlot);
    // dataPlotRef.zero();
    // ioMatrix::read("ImpactingBar.ref", "ascii", dataPlotRef);

    // double error = (dataPlot - dataPlotRef).normInf() ;
    // cout << "Error = " << error << endl;
    // if(error > 1e-11)
    // {
    //   std::cout << "Warning. The result is rather different from the reference file." << std::endl;
    //   std::cout << "Error = "<< (dataPlot - dataPlotRef).normInf()<<std::endl;
    //   return 1;
    // }


  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }


}
