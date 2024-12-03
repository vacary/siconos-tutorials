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
#include "myDS.h"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

user_defined::MyDS::MyDS(std::shared_ptr<siconos::algebra::SiconosVector> x0)
    : FirstOrderNonLinearDS(x0) {
  jacobianfVectorOver_x_ = std::make_shared<Matrix>(2, 2);
  _f = std::make_shared<Vector>(2);

  _M = std::make_shared<Matrix>(2, 2);
  _M->setZero();
  _M->setValue(0, 0, 1);
  _M->setValue(1, 1, 1);
}

void user_defined::MyDS::computefVector(
    const Eigen::Ref<siconos::algebra::SiconosVector> &state, double time) {
  // std::shared_ptr<siconos::algebra::SiconosVector> x=x();
  _f->setValue(0, -4.5 * x->getValue(0));
  _f->setValue(1, -1.5 * x->getValue(1));
  DEBUG_PRINT("MyDS::computeF");
  DEBUG_EXPR(x->display(););
  /*
  #ifdef SICONOS_DEBUG
    std::cout<<"MyDS::computeF with x="<<std::endl;
    x()->display();
    std::cout<<std::endl;
    std::cout<<"F(x)="<<std::endl;
    _f->display();
    std::cout<<std::endl;
  #endif
  */
}

void user_defined::MyDS::computeJacobianfOver_x(
    const Eigen::Ref<siconos::algebra::SiconosVector> &state, double time) {
  jacobianfVectorOver_x_->setValue(0, 0, -4.5);
  jacobianfVectorOver_x_->setValue(1, 0, 0);
  jacobianfVectorOver_x_->setValue(0, 1, 0);
  jacobianfVectorOver_x_->setValue(1, 1, -1.5);

  /*
  #ifdef SICONOS_DEBUG
    std::cout<<"MyDS::computeJacobianfx."<<std::endl;
  std::cout<<"Nabla f="<<std::endl;
    jacobianfVectorOver_x_->display();
    std::cout<<std::endl;
  #endif
  */
}

// void MyDS::computeRhs(double t)
// {
//   ;
// }
