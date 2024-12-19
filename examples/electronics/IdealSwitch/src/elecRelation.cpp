/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "elecRelation.h"

#include "circuit.h"

user_defined::elecRelation::elecRelation() : siconos::modeling::FirstOrderNonLinearR() {
  setComputehFunction([this](const siconos::algebra::BlockVector& state, double time,
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
                             Eigen::Ref<siconos::algebra::MapVectorType> y) {
    y.setZero();
#ifdef CLSC_CIRCUIT
    y(0) = lambda(4) - source(time);
    y(1) = state(0) - (lambda(3)) / sR;
    y(2) = lambda(2) - 20 + lambda(0) * (lambda(6) + sR1s);
    y(3) = lambda(2) + lambda(1) * (lambda(8) + sR1d);
    y(4) = state(0) - lambda(0) - lambda(1);
    y(5) = sR2 - lambda(6) - sR1s;
    y(6) = sAmpli * (lambda(4) - lambda(3)) + lambda(5);
    y(7) = sR2 - lambda(8) - sR1d;
    y(8) = -lambda(2) + lambda(7);
#else
    y(0) = -lambda(0) + source(t);
    y(1) = -lambda(0) + state(0) + (lambda(3) + sR1) * lambda(1);
    y(2) = sR2 - lambda(3) - sR1;
    y(3) = lambda(0) + lambda(2);
#endif
  });

  setComputegFunction([](const siconos::algebra::BlockVector& state, double time,
                         const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
                         siconos::algebra::BlockVector& result) {
#ifdef CLSC_CIRCUIT
    result(0) = (lambda(2) - lambda(3)) / sL;
#else
    result(0) = lambda(1) / sC;
#endif
  });

  setComputeJacobianhOver_stateFunction(
      [](const siconos::algebra::BlockVector& state, double time,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
         Eigen::Ref<siconos::algebra::MapType> result) {
        result.setZero();
#ifdef CLSC_CIRCUIT
        result(1, 0) = 1.;
        result(4, 0) = 1.;
#else
        result(1, 0) = 1.;
#endif
      });

  setComputeJacobianhOver_lambdaFunction(
      [](const siconos::algebra::BlockVector& state, double time,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
         Eigen::Ref<siconos::algebra::MapType> result) {
#ifdef CLSC_CIRCUIT
        result.setZero();

        result(0, 4) = 1;
        result(1, 3) = -1 / sR;
        result(2, 0) = lambda(6) + sR1s;
        result(2, 2) = 1;
        result(2, 6) = lambda(0);
        result(3, 1) = lambda(8) + sR1d;
        result(3, 2) = 1;
        result(3, 8) = lambda(1);

        result(4, 0) = -1;
        result(4, 1) = -1;
        result(5, 6) = -1;
        result(6, 3) = -sAmpli;
        result(6, 4) = sAmpli;
        result(6, 5) = 1;
        result(7, 8) = -1;
        result(8, 2) = -1;
        result(8, 7) = 1;
#else
        result(0, 0) = -1;
        result(1, 0) = -1;
        result(3, 0) = 1;
        result(5, 0) = lambda(3) + sR1;
        result(2, 1) = 1;
        result(4, 1) = lambda(1);
        result(5, 1) = -1;
#endif
      });

  setComputeJacobiangOver_lambdaFunction(
      [](const siconos::algebra::BlockVector& state, double time,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
         Eigen::Ref<siconos::algebra::MapType> result) {
        result.setZero();
#ifdef CLSC_CIRCUIT
        result(0, 2) = 1 / sL;
        result(0, 3) = -1 / sL;
#else
        result(0, 1) = 1 / sC;
#endif
      });
}

double user_defined::elecRelation::source(double t) {
  double daux = 0;
#ifdef CLSC_CIRCUIT
  double numT = t / sT;
  int aux = (int)floor(numT);
  daux = sE_plus - ((sE_plus - sE_moins) / sT) * t + (sE_plus - sE_moins) * aux;
  return daux;
#else
  daux = sin(sW * t);
  return daux;
#endif
}
