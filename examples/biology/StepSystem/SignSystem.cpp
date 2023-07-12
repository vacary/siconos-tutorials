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
//-----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <string>

#include "NonlinearRelationWithSign.hpp"
#include "NonlinearRelationWithSignInversed.hpp"
#include "SolverOptions.h"
#include "const.h"
#include "myDS.h"

using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char *argv[]) {
  try {
    // printf("argc %i\n", argc);
    int cmp = 0;
    int dimX = 2;

    char filename[50] = "sign.";

    //***** Set the initial condition
    // default, x0 = (1, 6)
    // else read from command line
    auto xti = std::make_shared<Vector>(dimX);
    if (argc == 1) {
      xti->setValue(0, 1);
      xti->setValue(1, 6);
      strncpy(&filename[5], "1.6.log", 7);
    } else if (argc == 3) {
      // printf("argv[0] %s\n", argv[0]);
      printf("xti(0) is set to %f\n", atof(argv[1]));
      printf("xti(1) is set to %f\n", atof(argv[2]));

      xti->setValue(0, atof(argv[1]));
      xti->setValue(1, atof(argv[2]));
      int sizeofargv1 = strlen(argv[1]);
      // printf("sizeofargv1 %i\n",sizeofargv1);
      strncpy(&filename[5], argv[1], sizeofargv1);
      int sizeofargv2 = strlen(argv[2]);
      // printf("sizeofargv2 %i\n",sizeofargv2);
      strncpy(&filename[5 + sizeofargv1], ".", 1);

      strncpy(&filename[5 + sizeofargv1 + 1], argv[2], sizeofargv2);
      strncpy(&filename[5 + sizeofargv1 + sizeofargv2 + 1], ".log", 4);

      printf("Output is written in filename %s\n", filename);
    } else {
      std::cout << "wrong  number of arguments = " << argc << "\n";
    }

    int NBStep = (int)floor(user_defined::sTf / user_defined::sStep);
    // NBStep =1;
    //*****BUILD THE DYNAMICAL SYSTEM

    auto aDS = std::make_shared<user_defined::MyDS>(xti);

    //******BUILD THE RELATION
    auto aR = std::make_shared<user_defined::NonlinearRelationWithSignInversed>();

    //*****BUILD THE NSLAW
    double ub = 1.;
    double lb = -1.;
    auto aNSL =
        std::make_shared<siconos::modeling::RelayNSL>(user_defined::sNSLawSize, lb, ub);

    //****BUILD THE INTERACTION
    auto aI = std::make_shared<siconos::modeling::Interaction>(aNSL, aR);

    //****BUILD THE model
    auto aN =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0, user_defined::sTf);
    aN->insertDynamicalSystem(aDS);
    aN->link(aI, aDS);

    // -- (1) OneStepIntegrators --
    auto aEulerMoreauOSI = std::make_shared<siconos::integrators::EulerMoreauOSI>(0.5);

    // -- (2) Time discretisation --
    auto aTD =
        std::make_shared<siconos::simulation::TimeDiscretisation>(0, user_defined::sStep);

    // -- (3) Non smooth problem
    auto osnspb =
        std::make_shared<siconos::nonsmooth_formulations::Relay>(SICONOS_RELAY_LEMKE);
    osnspb->numericsSolverOptions()->dparam[0] = 1e-08;

    // auto osnspb(new Relay(SICONOS_RELAY_ENUM);
    // osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS]=0;
    // // Multiple solutions 0 or 1 osnspb->numericsSolverOptions()->iparam[3]=48;

    osnspb->setNumericsVerboseMode(0);

    // -- (4) Simulation setup with (1) (2) (3)
    auto aS =
        std::make_shared<siconos::simulation::TimeStepping>(aN, aTD, aEulerMoreauOSI, osnspb);

    aS->setComputeResiduY(true);
    aS->setUseRelativeConvergenceCriteron(false);
    aS->setNewtonTolerance(5e-4);
    aS->setNewtonMaxIteration(20);
    aS->setResetAllLambda(false);
    // BUILD THE STEP INTEGRATOR

    auto x = aDS->x();
    auto vectorfield = aDS->rhs();
    auto y = aI->y(0);
    auto lambda = aI->lambda(0);

    unsigned int outputSize = 9;
    Matrix dataPlot(NBStep + 1, outputSize);

    std::cout << "=== Start of simulation: " << NBStep << " steps ===\n";

    printf("=== Start of simulation: %d steps ===  \n", NBStep);

    dataPlot(0, 0) = aN->t0();
    dataPlot(0, 1) = x->getValue(0);
    dataPlot(0, 2) = x->getValue(1);
    dataPlot(0, 3) = lambda->getValue(0);
    dataPlot(0, 4) = lambda->getValue(1);
    dataPlot(0, 5) = lambda->getValue(2);
    dataPlot(0, 6) = lambda->getValue(3);
    dataPlot(0, 7) = vectorfield->getValue(0);
    dataPlot(0, 8) = vectorfield->getValue(1);

    auto start = std::chrono::system_clock::now();
    for (int k = 0; k < NBStep; k++) {
#ifdef SICONOS_DEBUG
      std::std::cout << "-> Running step:" << k << "\n";
#endif
      cmp++;
      aS->advanceToEvent();

      dataPlot(cmp, 0) = aS->nextTime();
      dataPlot(cmp, 1) = x->getValue(0);
      dataPlot(cmp, 2) = x->getValue(1);
      dataPlot(cmp, 3) = lambda->getValue(0);
      dataPlot(cmp, 4) = lambda->getValue(1);
      dataPlot(cmp, 5) = lambda->getValue(2);
      dataPlot(cmp, 6) = lambda->getValue(3);

      aDS->computeRhs(aS->nextTime());

      if (cmp == 1)  // tricks just for display to avoid the computation of the initial Rhs
      {
        dataPlot(cmp - 1, 7) = vectorfield->getValue(0);
        dataPlot(cmp - 1, 8) = vectorfield->getValue(1);
      }

      dataPlot(cmp, 7) = vectorfield->getValue(0);
      dataPlot(cmp, 8) = vectorfield->getValue(1);

      aS->nextStep();

      // (*fout)<<cmp<<" "<<x->getValue(0)<<" "<<x->getValue(1)<<" "<<lambda->getValue(0)<<"
      // "<<lambda->getValue(1)<<" "<<lambda->getValue(2)<<" "<<lambda->getValue(3)<<"\n";
    }
    std::cout << "Computational time = "
              << "\n";
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";

    // dataPlot.resize(cmp, outputSize);
    siconos::algebra::io::write(filename, dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    if (argc == 1) {
      // Comparison with a reference file
      double error = 0.0, eps = 1e-11;
      if ((error = siconos::algebra::io::compareRefFile(dataPlot, "sign.1.6.ref", eps)) > eps)
        return 1;
    }

  } catch (...) {
    siconos::exception::process();
    return 1;
  }

  return 0;
}
