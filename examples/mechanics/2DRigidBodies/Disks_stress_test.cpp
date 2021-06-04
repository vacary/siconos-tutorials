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

/*!\file BouncingBallTS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/
#include <chrono>
#include <SiconosBodies.hpp>
#include <SiconosKernel.hpp>

#include <RigidBody2dDS.hpp>
#include "Contact2dR.hpp"

using namespace std;

static void compute_contact(SP::Contact2dR relation,
                            SP::RigidBody2dDS sphere_1,
                            SP::RigidBody2dDS sphere_2,
                            double R)
{
          SP::SiconosVector pc1 = relation->pc1();
          SP::SiconosVector pc2 = relation->pc2();
          SP::SiconosVector nc = relation->nc();

          SiconosVector& q_1 = *sphere_1->q();
          SiconosVector& q_2 = *sphere_2->q();

          double x_1 = q_1(0);
          double y_1 = q_1(1);
          double x_2 = q_2(0);
          double y_2 = q_2(1);

          double norm = sqrt((x_2-x_1)*(x_2-x_1)+(y_2-y_1)*(y_2-y_1)) ;
          (*nc)(0) = (x_2-x_1)/norm;
          (*nc)(1) = (y_2-y_1)/norm;

          (*pc1)(0) = x_1 + R* (*nc)(0) ;
          (*pc1)(1) = y_1 + R* (*nc)(1) ;
          (*pc2)(0) = x_2 - R* (*nc)(0) ;
          (*pc2)(1) = y_2 - R* (*nc)(1) ;


          //std::cout << "  distance  : "  << norm -2*R << std::endl;
          // std::cout << "  x1 x2 y1 y2  : "  << x_1 << " " << x_2 << " "<< y_1 << " " << y_2 << std::endl;
          // std::cout << "  normal  : ["  << (*nc)(0) << "," << (*nc)(1) <<"]" << std::endl;
          // std::cout << "  contact point  :"  << (*pc1)(0) << "," << (*pc1)(1) <<"]" << std::endl;
}
static void compute_contact_ground(SP::Contact2dR relation,
                                   SP::RigidBody2dDS sphere,
                                   double R)
{
          SP::SiconosVector pc1 = relation->pc1();
          SP::SiconosVector pc2 = relation->pc2();
          SP::SiconosVector nc = relation->nc();

          SiconosVector& q_1 = *sphere->q();

          double x_1 = q_1(0);
          double y_1 = q_1(1);

          double norm = sqrt((y_1)*(y_1)) ;
          (*nc)(0) = 0.0;
          (*nc)(1) = (y_1)/norm;
          (*pc1)(0) = x_1 - R* (*nc)(0) ;
          (*pc1)(1) = y_1 - R* (*nc)(1) ;
          (*pc2)(0) = x_1 ;
          (*pc2)(1) = R ;


          // std::cout << "  distance  : "  << norm << std::endl;
          // std::cout << "  x1 x2 y1 y2  : "  << x_1  << " "<< y_1 << std::endl;
          // std::cout << "  normal  : ["  << (*nc)(0) << "," << (*nc)(1) <<"]" << std::endl;
          // std::cout << "  contact point  :"  << (*pc1)(0) << "," << (*pc1)(1) <<"]" << std::endl;
          // std::cout << "  contact point  :"  << (*pc2)(0) << "," << (*pc2)(1) <<"]" << std::endl;
}
template<typename Duration = std::chrono::microseconds,
         typename F,
         typename ... Args>
typename Duration::rep profile(F&& fun,  Args&&... args)
{
  const auto beg = std::chrono::high_resolution_clock::now();
  //std::forward<F>(fun)(std::forward<Args>(args)...);
  std::invoke(fun, std::forward<Args>(args)...);
  const auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<Duration>(end - beg).count();
}

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time

    double h = 1e-2;// time step
    double T = 10*h;                  // final computation time
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    double R = 0.5; // Ball radius
    double m = 1; // Ball mass
    double g = 10; // Gravity

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." <<  endl;

    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 2. / 5 * m * R * R;

    double initialGap= 0.0;// - R/100.0;

    unsigned int n_col=500;
    unsigned int n_row=500;
    unsigned int n_spheres = n_col * n_row;

    std::vector<SP::SiconosVector> q0(n_spheres);
    std::vector<SP::SiconosVector> v0(n_spheres);

    for (unsigned int i =0 ; i < n_row; i++)
    {
       for (unsigned int j =0 ; j < n_col; j++)
       {
         (q0[j+ n_col*i]).reset(new SiconosVector(nDof));
         (v0[j+ n_col*i]).reset(new SiconosVector(nDof));
         (q0[j+ n_col*i])->setValue(0, position_init + i *2 * R + initialGap);
         (q0[j+ n_col*i])->setValue(1, position_init + j *2 * R + initialGap);

         (v0[j+ n_col*i])->setValue(0, velocity_init);
         (v0[j+ n_col*i])->setValue(1, velocity_init);

       }
    }

    SP::NonSmoothDynamicalSystem bouncingBall(new NonSmoothDynamicalSystem(t0, T));

    std::vector<SP::RigidBody2dDS> spheres(n_spheres);

    for (unsigned int i =0 ; i < n_row; i++)
    {
      for (unsigned int j =0 ; j < n_col; j++)
      {
        spheres[j+ n_col*i].reset(new RigidBody2dDS(q0[j+ n_col*i], v0[j+ n_col*i], Mass));
        SP::SiconosVector weight(new SiconosVector(nDof));
        (*weight)(1) = -10;//-m * g;
        // -- Set external forces (weight) --
        spheres[j+ n_col*i]->setFExtPtr(weight);
        bouncingBall->insertDynamicalSystem(spheres[j+n_col*i]);
      }
    }


    // --------------------
    // --- Interactions ---
    // --------------------
    // -- nslaw --
    double e = 0.9;

    // // Interaction ball-floor
    // //
    // SP::SimpleMatrix H(new SimpleMatrix(1, nDof));
    // (*H)(0, 0) = 1.0;
    // SP::Relation relation(new LagrangianLinearTIR(H));

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(e,0.0,0.1,2));



    std::vector<SP::Contact2dR > relation(n_spheres*(4));
    std::vector<SP::Interaction > interaction(n_spheres*(4));

    unsigned int id_interaction =0;
    for (unsigned int i =0 ; i < n_row-1; i++)
    {
      for (unsigned int j =0 ; j < n_col-1; j++)
      {
        unsigned int id_sphere= i + n_row*j;


        relation[id_interaction].reset(new Contact2dR());
        interaction[id_interaction].reset(new Interaction(nslaw, relation[id_interaction]));
        bouncingBall->link(interaction[id_interaction], spheres[id_sphere], spheres[id_sphere +1]);
        id_interaction++;
        relation[id_interaction].reset(new Contact2dR());
        interaction[id_interaction].reset(new Interaction(nslaw, relation[id_interaction]));
        bouncingBall->link(interaction[id_interaction], spheres[id_sphere], spheres[id_sphere + n_row]);
        id_interaction++;
        relation[id_interaction].reset(new Contact2dR());
        interaction[id_interaction].reset(new Interaction(nslaw, relation[id_interaction]));
        bouncingBall->link(interaction[id_interaction], spheres[id_sphere], spheres[id_sphere +1 + n_row]);
        id_interaction++;
      }
    }
    std::vector<SP::Contact2dR > relation_ground(n_row);
    std::vector<SP::Interaction > interaction_ground(n_row);

    for (unsigned int i =0 ; i < n_row; i++)
    {
      relation_ground[i].reset(new Contact2dR());
      interaction_ground[i].reset(new Interaction(nslaw, relation_ground[i]));
      bouncingBall->link(interaction_ground[i], spheres[i*n_row]);
    }


    // -------------
    // --- Model ---
    // -------------



    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));


    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::FrictionContact osnspb(new FrictionContact(2));
    osnspb->setMStorageType(NM_SPARSE);
    osnspb->setAssemblyType(REDUCED_DIRECT);

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(bouncingBall, t, OSI, osnspb));

    // =========================== End of model definition ===========================

    // ================================= Computation =================================


    int N = ceil((T - t0) / h)+1; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    SimpleMatrix dataPlot(N + 1, outputSize);

    // SP::SiconosVector q = ball->q();
    // SP::SiconosVector v = ball->velocity();
    // SP::SiconosVector p = ball->p(1);
    // SP::SiconosVector lambda = inter->lambda(1);
    // SP::SiconosVector q1 = ball1->q();
    // SP::SiconosVector v1 = ball->velocity();
    // SP::SiconosVector p1 = ball->p(1);
    // SP::SiconosVector lambda1 = inter1->lambda(1);


    // dataPlot(0, 0) = bouncingBall->t0();
    // dataPlot(0, 1) = (*q)(0);
    // dataPlot(0, 2) = (*v)(0);
    // dataPlot(0, 3) = (*p)(0);
    // dataPlot(0, 4) = (*lambda)(0);
    // dataPlot(0, 5) = (*q1)(0);
    // dataPlot(0, 6) = (*v1)(0);
    // dataPlot(0, 7) = (*p1)(0);
    // dataPlot(0, 8) = (*lambda1)(0);
    // // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    std::chrono::time_point<std::chrono::system_clock> start, end, before, after;
    start = std::chrono::system_clock::now();
    int elapsed_nextStep =0,  elapsed_computeOneStep=0;
    int
      time_initialize =0,
      time_initializeOSIAssociations =0,
      time_initializeIndexSets=0,
      time_updateWorldFromDS =0,
      time_initializeNSDSChangelog =0,
      time_updateOutput=0,
      time_initOSNS=0,
      time_firstInitialize =0,
      time_initializeNewtonLoop =0,
      time_prepareNewtonIteration=0,
      time_computeFreeState =0,
      time_computeOneStepNSProblem=0,
      time_precompute =0,
      time_updateMu =0,
      time_solve=0,
      time_postCompute=0,
      time_updateInput=0,
      time_updateState=0;
    // SP::SiconosVector rpc = relation->relPc1();
    // SP::SiconosVector nnc = relation->relNc();

    // SP::SiconosVector rpc1 = relation1->relPc1();
    // SP::SiconosVector rpc2 = relation1->relPc2();
    // SP::SiconosVector nnc1 = relation1->relNc();

    const Simulation & simulation = *s;


    while(s->hasNextEvent())
    {


      std::cout << "\n\n\nIteration :" << k << std::endl;
      osnspb->setNumericsVerboseMode(true);
      id_interaction =0;
      for (unsigned int i =0 ; i < n_row-1; i++)
      {
        for (unsigned int j =0 ; j < n_col-1; j++)
        {
          unsigned int id_sphere= i + n_row*j;
          //std::cout << "interaction id :"  << id_interaction  <<  std::endl;
          compute_contact(relation[id_interaction], spheres[id_sphere], spheres[id_sphere+1],R);
          id_interaction++;
          //std::cout << "interaction id :"  << id_interaction  <<  std::endl;
          compute_contact(relation[id_interaction], spheres[id_sphere], spheres[id_sphere+n_row],R);
          id_interaction++;
          //std::cout << "interaction id :"  << id_interaction  <<  std::endl;
          compute_contact(relation[id_interaction], spheres[id_sphere], spheres[id_sphere+1 + n_row],R);
          id_interaction++;
        }
      }
      for (unsigned int i =0 ; i < n_row; i++)
      {
        compute_contact_ground(relation_ground[i], spheres[i*n_row], R);
      }

      //getchar();


      // // a fake contact detection
      // (*rpc)(0) = -R;
      // (*rpc)(1) = 0.0;
      // (*nnc)(0) = 1.0;
      // (*nnc)(1) = 0.0;

      // (*rpc1)(0) = -R;
      // (*rpc1)(1) = 0.0;

      // (*rpc2)(0) = R;
      // (*rpc2)(1) = 0.0;
      // (*nnc1)(0) = 1.0;
      // (*nnc1)(1) = 0.0;


      int time=0;
      osnspb->setNumericsVerboseMode(true);
      // time = profile<std::chrono::microseconds>(&TimeStepping::computeOneStep,s);
      // std::cout << "time elapsed in computeOneStep" << time<< std::endl;
      // elapsed_computeOneStep +=  time;

      //time  =profile<std::chrono::microseconds>(&TimeStepping::initialize,s);
      //time_initialize+= time;
      time  =profile<std::chrono::microseconds>(&Simulation::initializeOSIAssociations,s);
      time_initializeOSIAssociations += time;
      time  =profile<std::chrono::microseconds>(&Simulation::initializeIndexSets,s);
      time_initializeIndexSets += time;
      time  =profile<std::chrono::microseconds>(&Simulation::updateWorldFromDS,s);
      time_updateWorldFromDS += time;
      time  =profile<std::chrono::microseconds>(&Simulation::initializeNSDSChangelog,s);
      time_initializeNSDSChangelog += time;

      time  =profile<std::chrono::microseconds>(&Simulation::updateOutput,s,0);
      time_updateOutput += time;
      time  =profile<std::chrono::microseconds>(&Simulation::initOSNS,s);
      time_initOSNS += time;
      time  =profile<std::chrono::microseconds>(&Simulation::firstInitialize,s);
      time_firstInitialize += time;


      time  =profile<std::chrono::microseconds>(&TimeStepping::initializeNewtonLoop,s);
      time_initializeNewtonLoop += time;
      time  =profile<std::chrono::microseconds>(&TimeStepping::prepareNewtonIteration,s);
      time_prepareNewtonIteration += time;
      time  =profile<std::chrono::microseconds>(&TimeStepping::computeFreeState,s);
      time_computeFreeState += time;
      // time  =profile<std::chrono::microseconds>(&TimeStepping::computeOneStepNSProblem,s,SICONOS_OSNSP_TS_VELOCITY);
      // time_computeOneStepNSProblem += time;
      time  =profile<std::chrono::microseconds>(&FrictionContact::preCompute,osnspb,0.0);
      time_precompute += time;
      time  =profile<std::chrono::microseconds>(&FrictionContact::updateMu,osnspb);
      time_updateMu += time;
      time  =profile<std::chrono::microseconds>(&FrictionContact::solve,osnspb,SP::FrictionContactProblem());
      time_solve += time;
      time  =profile<std::chrono::microseconds>(&FrictionContact::postCompute,osnspb);
      time_postCompute += time;

      time  =profile<std::chrono::microseconds>(&TimeStepping::updateInput,s,0);
      time_updateInput += time;
      time  =profile<std::chrono::microseconds>(&TimeStepping::updateState,s,0);
      time_updateState += time;


      std::cout<< "size output " << osnspb->getSizeOutput() << std::endl;;
      // // --- Get values to be plotted ---
      // dataPlot(k, 0) =  s->nextTime();
      // dataPlot(k, 1) = (*q)(0);
      // dataPlot(k, 2) = (*v)(0);
      // dataPlot(k, 3) = (*p)(0);
      // dataPlot(k, 4) = (*lambda)(0);
      // dataPlot(k, 5) = (*q1)(0);
      // dataPlot(k, 6) = (*v1)(0);
      // dataPlot(k, 7) = (*p1)(0);
      // dataPlot(k, 8) = (*lambda1)(0);



      elapsed_nextStep += profile<std::chrono::microseconds>(&TimeStepping::nextStep,s);
      k++;

    }
    cout  << "End of computation - Number of iterations done: " << k - 1 << endl;
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    cout << "Computation time : " << elapsed << " ms" << endl;

    cout << "time Computeonestep          : " <<  elapsed_computeOneStep/1000 << " ms" << endl<< endl;;

    cout << "time initialize              : " <<  time_initialize/1000 << " ms" << endl;
    cout << "   time_initializeOSIAssociations : " << time_initializeOSIAssociations/1000 << " ms" << endl;
    cout << "   time_initializeIndexSets       : "<<time_initializeIndexSets/1000 << " ms" << endl;
    cout << "   time_updateWorldFromDS         : "<< time_updateWorldFromDS/1000 << " ms" << endl;
    cout << "   time_initializeNSDSChangelog   : "<< time_initializeNSDSChangelog/1000 << " ms" << endl;
    cout << "   time_updateOutput              : "<< time_updateOutput/1000 << " ms" << endl;
    cout << "   time_initOSNS                  : "<< time_initOSNS/1000 << " ms" << endl;
    cout << "   time_firstInitialize           : " << time_firstInitialize/1000 << " ms" << endl;


    
    cout << "time initializeNewtonLoop    : " <<  time_initializeNewtonLoop/1000 << " ms" << std::endl;
    cout << "time prepareNewtonIteration  : " <<  time_prepareNewtonIteration/1000 << " ms" << std::endl;
    cout << "time computeFreeState        : " <<  time_computeFreeState/1000  << " ms" << std::endl;
    cout << "time computeOneStepNSProblem : " <<  time_computeOneStepNSProblem/1000 << " ms" << std::endl;
    cout << "    time precompute               : " <<   time_precompute/1000 << " ms" << std::endl;
    cout << "    time updateMu                 : " <<  time_updateMu/1000 << " ms" << std::endl;
    cout << "    time solve                    : " <<  time_solve/1000 << " ms" << std::endl;
    cout << "    time postCompute              : " <<  time_postCompute/1000 << " ms" << std::endl;


    cout << "time updateInput             : " <<  time_updateInput/1000 << " ms" << std::endl;
    cout << "time updateState             : " <<  time_updateState/1000  << " ms" << std::endl << endl;
    cout << "time nextStep                : " <<  elapsed_nextStep/1000 << " ms" << endl;


  //   // --- Output files ---
  //   cout << "====> Output file writing ..." << endl;
  //   dataPlot.resize(k, outputSize);
  //   ioMatrix::write("Ball2D.dat", "ascii", dataPlot, "noDim");
  //   double error=0.0, eps=1e-12;
  //   if((error=ioMatrix::compareRefFile(dataPlot, "Ball2D.ref", eps)) >= 0.0
  //       && error > eps)
  //     return 1;

  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }



}
