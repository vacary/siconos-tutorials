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
    DEBUG_EXPR(x.display());
    DEBUG_EXPR(lambda.display());
    y.setValue(0, 4.0 - x(0));
    y.setValue(1, 4.0 - x(1));
    y.setValue(2, 8.0 - x(0));
    y.setValue(3, 8.0 - x(1));
    DEBUG_EXPR(y.display());
  });

  setComputegFunction([](const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
                         siconos::algebra::BlockVector& res) {
    DEBUG_PRINTF("user_defined::NonlinearRelation::computeg at time %e\n ", t);
    DEBUG_EXPR(lambda.display());

    res.setValue(0, 40.0 * (1 - lambda(2)) * (lambda(1)));
    res.setValue(1, 40.0 * (lambda(0)) * (1 - lambda(3)));
    /*
    #ifdef SICONOS_DEBUG
      std::cout<<"user_defined::NonlinearRelation::computeg with lambda="<<std::endl;
      lambda.display();
      std::cout<<std::endl;
      std::cout<<"user_defined::NonlinearRelation::computeg modif g_alpha : \n";
      inter.data(g_alpha)->display();
      std::cout<<std::endl;
    #endif
    */
    DEBUG_EXPR(res.display());
  });

  setComputeJacobianhOver_stateFunction(
      [](const siconos::algebra::BlockVector& state,
         const Eigen::Ref<const siconos::algebra::SiconosVector>& lam) {
        DEBUG_PRINTF("user_defined::NonlinearRelation:: compute jacobian h over lambda  at time %e\n ", t);

        jacobianhOver_state_view_.setValue(0, 0, -1);
        jacobianhOver_state_view_.setValue(0, 1, 0);
        jacobianhOver_state_view_.setValue(1, 0, 0);
        jacobianhOver_state_view_.setValue(1, 1, -1);
        jacobianhOver_state_view_.setValue(2, 0, -1);
        jacobianhOver_state_view_.setValue(2, 1, 0);
        jacobianhOver_state_view_.setValue(3, 0, 0);
        jacobianhOver_state_view_.setValue(3, 1, -1);
        DEBUG_EXPR(std::cout << jacobianhOver_state_view_ << "\n";);
      });

  setComputeJacobiangOver_lambdaFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
        DEBUG_PRINTF("user_defined::NonlinearRelation::compute jacobian g over lambda at time %e\n ", t);
        DEBUG_EXPR(lambda.display());

        jacobiangOver_lambda_view_.setValue(0, 0, 0);
        jacobiangOver_lambda_view_.setValue(1, 0, 40.0 * (1 - lambda(3)));

        jacobiangOver_lambda_view_.setValue(0, 1, 40.0 * (1 - lambda(2)));
        jacobiangOver_lambda_view_.setValue(1, 1, 0);

        jacobiangOver_lambda_view_.setValue(0, 2, -40.0 * lambda(1));
        jacobiangOver_lambda_view_.setValue(1, 2, 0);

        jacobiangOver_lambda_view_.setValue(0, 3, 0);
        jacobiangOver_lambda_view_.setValue(1, 3, -40.0 * lambda(0));
        DEBUG_EXPR(jacobiangOver_lambda_view_.display());
      });
}
