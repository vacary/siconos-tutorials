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
#include "NonlinearRelation.hpp"

// #include "const.h"

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

user_defined::NonlinearRelation::NonlinearRelation() : FirstOrderType2R{} {
  setComputehFunction([](const siconos::algebra::BlockVector& state,
                         const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
                         Eigen::Ref<siconos::algebra::SiconosVector> y) {
    DEBUG_PRINTF("user_defined::NonlinearRelation::computeh at time %e\n ", t);
    DEBUG_EXPR(siconos::algebra::print(x));
    DEBUG_EXPR(siconos::algebra::print(lambda));
    y(0) = 4.0 - state(0);
    y(1) = 4.0 - state(1);
    y(2) = 8.0 - state(0);
    y(3) = 8.0 - state(1);
    DEBUG_EXPR(siconos::algebra::print(y));
  });

  setComputegFunction([](const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
                         siconos::algebra::BlockVector& res) {
    DEBUG_EXPR(siconos::algebra::print(lambda));

    res(0) = 40.0 * (1 - lambda(2)) * (lambda(1));
    res(1) = 40.0 * (lambda(0)) * (1 - lambda(3));

    DEBUG_EXPR(siconos::algebra::print(res));
  });

  setComputeJacobianhOver_stateFunction(
      [](const siconos::algebra::BlockVector& state,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& lam,
         Eigen::Ref<siconos::algebra::MapType> result) {
        result.setZero();
        result.setValue(0, 0, -1);
        result.setValue(1, 1, -1);
        result.setValue(2, 0, -1);
        result.setValue(3, 1, -1);
        DEBUG_EXPR(std::cout << result << "\n";);
      });

  setComputeJacobiangOver_lambdaFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
         Eigen::Ref<siconos::algebra::MapType> result) {
        DEBUG_PRINTF(
            "user_defined::NonlinearRelation::compute jacobian g over lambda at time %e\n ",
            t);
        DEBUG_EXPR(siconos::algebra::print(lambda));
        result.setValue(1, 0, 40.0 * (1 - lambda(3)));
        result.setValue(0, 1, 40.0 * (1 - lambda(2)));
        result.setValue(0, 2, -40.0 * lambda(1));
        result.setValue(1, 3, -40.0 * lambda(0));
        DEBUG_EXPR(siconos::algebra::print(result));
      });
}
