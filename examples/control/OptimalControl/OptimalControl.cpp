/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include <NumericsVerbose.h>
#include <stdio.h>
#include <stdlib.h>

#include <SiconosKernel.hpp>
#include <adjointInput.hpp>
#include <chrono>

#include "myDS.h"

using namespace std;
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

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
  int cmp = 0;

  /************************************************************/
  /************************************************************/

  int dimX = 4;

  auto x0 = std::make_shared<Vector>(dimX);

  // Point de départ hors arc singulier
  (*x0)(0) = 3.4999939172;
  (*x0)(1) = -2.2788416237;
  (*x0)(2) = 1.1935988302;
  (*x0)(3) = -0.6365413023;

  double sT = 10;
  double sStep = 2e-3;
  unsigned int NBStep = floor(sT / sStep);
  // NBStep =2;

  //  NBStep = 3;
  //*****BUILD THE DYNAMIC SYSTEM
  auto aDS = std::make_shared<user_defined::MyDS>(*x0);

  //******BUILD THE RELATION
  auto aR = std::make_shared<user_defined::adjointInput>();

  int sN = 2;

  //*****BUILD THE NSLAW
  auto aNSL = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(sN);

  //****BUILD THE INTERACTION
  auto aI = std::make_shared<siconos::modeling::Interaction>(aNSL, aR);
  //****BUILD THE SYSTEM
  auto aM = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0, sT);
  aM->insertDynamicalSystem(aDS);
  aM->link(aI, aDS);
  auto aTD = std::make_shared<siconos::simulation::TimeDiscretisation>(0, sStep);
  auto aS = std::make_shared<siconos::simulation::TimeStepping>(aM, aTD);
  aS->setComputeResiduY(true);
  aS->setComputeResiduR(true);
  aS->setUseRelativeConvergenceCriteron(false);
  aS->setNewtonTolerance(1.1e-11);
  aS->setNewtonMaxIteration(50);
  //*****BUILD THE STEP INTEGRATOR
  auto aEulerMoreauOSI = std::make_shared<siconos::integrators::EulerMoreauOSI>(0.5);
  aS->insertIntegrator(aEulerMoreauOSI);

  //**** BUILD THE STEP NS PROBLEM
  auto aLCP = std::make_shared<siconos::nonsmooth_formulations::LCP>(SICONOS_LCP_ENUM);
  //  aLCP.reset=
  //  std::make_shared<siconos::nonsmooth_formulations::LCP>(SICONOS_LCP_NEWTONFB);

  aS->insertNonSmoothProblem(aLCP);

  numerics_set_verbose(0);

  auto x = aDS->x();
  auto y = aI->y(0);
  auto lambda = aI->lambda(0);

  unsigned int outputSize = 9;  // number of required data
  Matrix dataPlot(NBStep + 1, outputSize);

  auto z = aDS->x();

  dataPlot(0, 0) = aM->t0();  // Initial time of the model
  dataPlot(0, 1) = (*z)(0);
  dataPlot(0, 2) = (*z)(1);
  dataPlot(0, 3) = (*z)(2);
  dataPlot(0, 4) = (*z)(3);
  dataPlot(0, 5) = (*lambda)(0);
  dataPlot(0, 6) = (*lambda)(1);
  dataPlot(0, 7) = (*y)(0);
  dataPlot(0, 8) = (*y)(1);

  // do simulation while events remains in the "future events" list of events manager.
  cout << " ==== Start of  simulation : " << NBStep << " steps====\n";
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  unsigned int k = 0;
  while (aS->hasNextEvent()) {
    k++;
    //      if (cmp==150)
    // numerics_set_verbose(à);
    //      else if (cmp==151)
    numerics_set_verbose(0);

    cmp++;

    // solve ...
    aS->computeOneStep();
    x = aDS->x();
    lambda = aI->lambda(0);
    dataPlot(k, 0) = aS->nextTime();  // Initial time of the model
    dataPlot(k, 1) = (*x)(0);
    dataPlot(k, 2) = (*x)(1);
    dataPlot(k, 3) = (*x)(2);
    dataPlot(k, 4) = (*x)(3);
    dataPlot(k, 5) = (*lambda)(0);
    dataPlot(k, 6) = (*lambda)(1);
    dataPlot(k, 7) = (*y)(0);
    dataPlot(k, 8) = (*y)(1);
    aS->nextStep();
  }

  cout << "===== End of simulation. ==== \n";

  // --- Output files ---
  cout << "====> Output file writing ...\n";
  siconos::algebra::io::write("OptimalControl.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);

  double error = 0.0, eps = 1e-08;
  if ((error = siconos::algebra::io::compareRefFile(dataPlot, "OptimalControl.ref", eps)) >
      eps)
    return 1;

  return 0;
}
