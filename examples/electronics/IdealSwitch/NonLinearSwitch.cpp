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

#include <MLCP_Solvers.h>
#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <string>

#include "circuit.h"
#include "elecRelation.h"
#include "myDS.h"

using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

/************************************************************/
/************************************************************/
/************************************************************/
/*call back for the source*/
/*call back for the formulation with inversion*/
// void (bLDS) (double t, unsigned int N, double* b, unsigned int z, double*zz){
// }

/************************************************************/
/************************************************************/
/************************************************************/
/************************************************************/
/*main program*/

int main() {
  string solverDirEnum = "DIRECT_ENUM";
  string solverDirPath = "DIRECT_PATH";
  string solverDirSimplex = "DIRECT_SIMPLEX";
  string solverEnum = "ENUM";
  string solverSimplex = "SIMPLEX";
  string solverPath = "PATH";
  string* solverName = 0;
  bool diodeIsOn = true;
  bool switchIsOn = true;
  bool stateChanged = true;
  // One Step non smooth problem

  double* floatWorkingMem = 0;
  int* intWorkingMem = 0;

  // int freq = 1000;
  // int Nfreq = 0;
  int cmp = 0;

  // int NbDataMax = 10000;
  // int NData = 0;

  /************************************************************/
  /************************************************************/
  /*Solver options*/
  solverName = &solverEnum;

  int dimX = 1;
  // SimpleMatrix * M = 0;
  // SimpleMatrix * A = 0;
  // SiconosVector* As = 0;
  // SiconosVector* mti = 0;

  auto xti = std::make_shared<Vector>(dimX);
  xti->setValue(0, 0);

  int NBStep = (int)floor(user_defined::sTf / user_defined::sStep);

  // NBStep = 130;
  //*****BUILD THE DYNAMIC SYSTEM
  auto aDS = std::make_shared<user_defined::MyDS>(xti);

  //******BUILD THE RELATION
  // SimpleMatrix* C = 0;
  // SimpleMatrix* D = 0;
  // SimpleMatrix* B = 0;
  auto aR = std::make_shared<user_defined::elecRelation>();

  //*****BUILD THE NSLAW
  auto aNSL = std::make_shared<siconos::modeling::MixedComplementarityConditionNSL>(
      user_defined::sN, user_defined::sM);
  /*
    if (ACE_FORMULATION==ACE_FORMULATION_SEMI_EXPLICT){
    NSLawSize=m+s;
    aNSL.reset= std::make_shared<siconos::modeling::MixedComplementarityConditionNSL>(m,s);
    }else{
    NSLawSize=m;
    aNSL.reset= std::make_shared<siconos::modeling::ComplementarityConditionNSL>(m);
    }
  */

  //****BUILD THE INTERACTION
  auto aI = std::make_shared<siconos::modeling::Interaction>(aNSL, aR);
  //  aI->insert(LSDiodeBridge);
  //****BUILD THE SYSTEM

  auto aN =
      std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0, user_defined::sTf);
  aN->insertDynamicalSystem(aDS);
  aN->link(aI, aDS);

  // -- (1) OneStepIntegrators --
  auto aEulerMoreauOSI = std::make_shared<siconos::integrators::EulerMoreauOSI>(0.5);

  // -- (2) Time discretisation --
  auto aTD = std::make_shared<siconos::simulation::TimeDiscretisation>(0, user_defined::sStep);

  // -- (3) Non smooth problem
  auto aMLCP = std::make_shared<siconos::nonsmooth_formulations::MLCP>();

  // -- (4) Simulation setup with (1) (2) (3)
  auto aS =
      std::make_shared<siconos::simulation::TimeStepping>(aN, aTD, aEulerMoreauOSI, aMLCP);
  aS->setComputeResiduY(true);
  aS->setComputeResiduR(true);
  aS->setUseRelativeConvergenceCriteron(false);
  aS->setNewtonMaxIteration(20);
  aS->setNewtonTolerance(1e-11);
  aS->setResetAllLambda(false);
  // To compute necessary information for memory allocator

  aS->initialize();

  aMLCP->preCompute(0.0);
  cout << "nonSmoothDynamicalSystem()->isLinear() : " << boolalpha << aN->isLinear() << "\n";

  cout << "nonSmoothDynamicalSystem()->topology()->hasChanged() : " << boolalpha
       << aN->topology()->hasChanged() << "\n";

  //*****BUILD THE STEP INTEGRATOR
  //  SP::NonSmoothSolver  mySolver( new
  //  NonSmoothSolver((*solverName),iparam,dparam,floatWorkingMem,intWorkingMem);

  //**** BUILD THE STEP NS PROBLEM

  //      numerics_set_verbose(1);

  auto x = aDS->x();
  auto y = aI->y(0);
  auto lambda = aI->lambda(0);
  ofstream* fout = new ofstream("simu.log");
  fout->precision(10);
  ifstream* fin = new ifstream("NonLinearSwitch.ref");
  fin->precision(10);
  // unsigned int count = 0; // events counter.
  // do simulation while events remains in the "future events" list of events manager.
  cout << " ==== Start of  simulation : " << NBStep << " steps====\n";
#ifdef CLSC_CIRCUIT
#else
  (*fout) << "C_charge "
          << "V1 "
          << "R(t)\n";
#endif

  for (int k = 0; k < NBStep; k++) {
    //      if (cmp==150)
    //        numerics_set_verbose(1);
    //      else if (cmp==151)
    //        numerics_set_verbose(0);
    // cout << "..." << cmp << "\n";
    cmp++;
    // solve ...
    aS->computeOneStep();
    // aMLCP->display();
    aS->nextStep();
    x = aDS->x();
    lambda = aI->lambda(0);
#ifdef CLSC_CIRCUIT

    // std::cout<<"x="<<x->getValue(0)<<" Is="<<lambda->getValue(0)<<"
    // Id="<<lambda->getValue(1)<<" V3="<<lambda->getValue(2); std::cout<<"
    // V4="<<lambda->getValue(3)<<" V5="<<lambda->getValue(4)<<" l6="<<lambda->getValue(5)<<"
    // l7="<<lambda->getValue(6); std::cout<<" l8="<<lambda->getValue(7)<<"
    // l9="<<lambda->getValue(8)<<std::"\n";
    stateChanged = false;
    if (lambda->getValue(6) > 1) {
      if (switchIsOn || k == 0) {
        switchIsOn = false;
        stateChanged = true;
      }
    } else {
      if (!switchIsOn || k == 0) {
        switchIsOn = true;
        stateChanged = true;
      }
    }
    if (lambda->getValue(8) > 1) {
      if (diodeIsOn || k == 0) {
        diodeIsOn = false;
        stateChanged = true;
      }
    } else {
      if (!diodeIsOn || k == 0) {
        diodeIsOn = true;
        stateChanged = true;
      }
    }
    if (stateChanged) {
      if (switchIsOn)
        std::cout << "SWITCH=ON";
      else
        std::cout << "SWITCH=OFF";
      if (diodeIsOn)
        std::cout << " DIODE=ON";
      else
        std::cout << " DIODE=OFF";
      std::cout << "\n";
    }

    (*fout) << cmp << " " << x->getValue(0) << " " << lambda->getValue(0) << " "
            << lambda->getValue(1) << " " << lambda->getValue(2) << " " << lambda->getValue(3)
            << " " << lambda->getValue(4) << " " << lambda->getValue(5) << "\n";
    int cmpR;
    double xR;
    string sz;

    (*fin) >> cmpR >> xR;

    getline(*fin, sz);
    // cout << "==== difference = " <<fabs(xR - x->getValue(0))  <<"\n";

    if (fabs(xR - x->getValue(0)) > 10e-7) {
      cout << "==== simulation is stopped because of a too large difference with a referenced "
              "trajectory. ==== \n";
      cout << "==== difference = " << fabs(xR - x->getValue(0)) << "\n";
      // return 1;
    }
#else
    (*fout) << cmp << " " << x->getValue(0) << " " << lambda->getValue(0) << " "
            << lambda->getValue(3) + sR1 << "\n";
#endif
  }
  delete fout;
  delete fin;
  cout << "===== End of simulation. ==== \n";
  return 0;
}
