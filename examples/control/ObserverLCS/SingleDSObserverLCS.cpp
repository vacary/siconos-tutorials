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

#include <SiconosKernel.hpp>
#include <SiconosPointers.hpp>
#include <chrono>
#include <string>

using namespace std;
using SimpleMatrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;
using namespace std;

int main(int argc, char* argv[]) {
  // Exception handling
  try {
    // == User-defined parameters ==
    unsigned int ndof = 4;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 25;      // Total simulation time
    double h = 1.0e-3;  // Time step
    double Vinit = 10.0;
    unsigned int noutput = 1;

    // ================= Creation of the model =======================

    // == Creation of the NonSmoothDynamicalSystem ==
    // DynamicalSystem(s)
    SimpleMatrix A(2, 2);  // All components of A are automatically set to 0.
    A(0, 0) = 1.0;
    A(0, 1) = 1.0;
    A(1, 0) = 3.0;
    A(1, 1) = 1.0;
    A = 0.1 * A;
    SimpleMatrix TildeA(ndof, ndof);  // All components of A are automatically set to 0.
    TildeA(0, 0) = A(0, 0);
    TildeA(0, 1) = A(0, 1);
    TildeA(1, 0) = A(1, 0);
    TildeA(1, 1) = A(1, 1);

    SimpleMatrix L(2, noutput);
    L(0, 0) = 1.0;
    L(1, 0) = 1.0;
    L = 0.1 * L;
    SimpleMatrix G(noutput, 2);
    G(0, 0) = 2.0;
    G(0, 1) = 2.0;

    SimpleMatrix hatA(2, 2);
    hatA = A - prod(L, G);
    TildeA(2, 2) = hatA(0, 0);
    TildeA(2, 3) = hatA(0, 1);
    TildeA(3, 2) = hatA(1, 0);
    TildeA(3, 3) = hatA(1, 1);

    SimpleMatrix LG(2, 2);
    LG = prod(L, G);
    TildeA(2, 0) = LG(0, 0);
    TildeA(3, 0) = LG(1, 0);
    TildeA(2, 1) = LG(0, 1);
    TildeA(3, 1) = LG(1, 1);

    auto x0 = std::make_shared<Vector>(ndof);
    (*x0)(0) = Vinit;
    auto processObserver = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        x0, siconos::pointers::createSPtr(TildeA));
    processObserver->setComputebFunction("SingleDSObserverLCSPlugin", "computeU");

    // Relations
    unsigned int ninter = 2;  // dimension of your Interaction = size of y and lambda vectors
    SimpleMatrix B(ndof, ninter);
    B(0, 0) = -1.0;
    B(1, 0) = 1.0;
    B(2, 1) = -1.0;
    B(3, 1) = 1.0;
    SimpleMatrix C(ninter, ndof);
    C(0, 0) = -1.0;
    C(0, 1) = 1.0;
    C(1, 2) = -1.0;
    C(1, 3) = 1.0;

    auto myProcessRelation = std::make_shared<siconos::modeling::FirstOrderLinearR>(
        siconos::pointers::createSPtr(C),
        siconos::pointers::createSPtr(B));

    myProcessRelation->setComputeEFunction("SingleDSObserverLCSPlugin", "computeE");

    SimpleMatrix D(ninter, ninter);
    D(0, 0) = 1.0;
    D(1, 1) = 1.0;
    // myProcessRelation->setD(D);
    // return 0;

    // NonSmoothLaw
    unsigned int nslawSize = 2;
    auto myNslaw = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(nslawSize);

    auto myProcessInteraction =
        std::make_shared<siconos::modeling::Interaction>(myNslaw, myProcessRelation);

    // Model
    auto ObserverLCS = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    ObserverLCS->insertDynamicalSystem(processObserver);
    ObserverLCS->link(myProcessInteraction, processObserver);
    // TimeDiscretisation
    auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // == Creation of the Simulation ==
    auto s = std::make_shared<siconos::simulation::TimeStepping>(ObserverLCS, td);

    // OneStepIntegrator
    double theta = 0.5;
    // One Step Integrator
    auto myIntegrator = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);
    s->insertIntegrator(myIntegrator);

    // One Step non smooth problem

    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    s->insertNonSmoothProblem(osnspb);

    // ================================= Computation =================================

    int k = 0;                                // Current step
    unsigned int N = ceil((T - t0) / h) + 1;  // Number of time steps
    unsigned int outputSize = 10;             // number of required data
    SimpleMatrix dataPlot(N, outputSize);
    auto processLambda = myProcessInteraction->lambda(0);

    myProcessInteraction->computeOutput(t0, 0);
    // We get values for the initial time step:
    // time
    dataPlot(k, 0) = s->nextTime();
    ;
    dataPlot(k, 1) = (*processObserver->x())(0);  // Observer x(1)
    dataPlot(k, 2) = (*processObserver->x())(1);  // Observer x(2)
    dataPlot(k, 3) = (*processObserver->x())(2);  // Process x(1)
    dataPlot(k, 4) = (*processObserver->x())(3);  // Process x(2)
    dataPlot(k, 5) = (*processLambda)(0);
    dataPlot(k, 6) = (*processLambda)(1);
    dataPlot(k, 7) = (*processObserver->b())(0);
    dataPlot(k, 8) = abs((*processObserver->x())(0) - (*processObserver->x())(2));
    dataPlot(k, 9) = abs((*processObserver->x())(1) - (*processObserver->x())(3));

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    // Simulation loop
    while (s->hasNextEvent()) {
      k++;

      // get current time step

      s->computeOneStep();

      dataPlot(k, 0) = s->nextTime();
      ;
      dataPlot(k, 1) = (*processObserver->x())(0);
      dataPlot(k, 2) = (*processObserver->x())(1);
      dataPlot(k, 3) = (*processObserver->x())(2);
      dataPlot(k, 4) = (*processObserver->x())(3);
      dataPlot(k, 5) = (*processLambda)(0);
      dataPlot(k, 6) = (*processLambda)(1);
      dataPlot(k, 7) = (*processObserver->b())(0);
      dataPlot(k, 8) = abs((*processObserver->x())(0) - (*processObserver->x())(2));
      dataPlot(k, 9) = abs((*processObserver->x())(1) - (*processObserver->x())(3));

      s->nextStep();
    }

    // Write the results into the file "ObserverLCS.dat"
    siconos::algebra::io::write("SingleDSObserverLCS.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-9;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "SingleDSObserverLCS.ref",
                                                      eps)) > eps)
      return 1;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
