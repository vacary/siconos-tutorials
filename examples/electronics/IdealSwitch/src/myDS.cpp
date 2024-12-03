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

user_defined::MyDS::MyDS(std::shared_ptr<siconos::algebra::SiconosVector> x0)
    : FirstOrderNonLinearDS(x0) {
  jacobianfVectorOver_x_ = std::make_shared<siconos::algebra::SiconosMatrix>(1, 1);
  _f = std::make_shared<siconos::algebra::SiconosVector>(1);
  _M = std::make_shared<siconos::algebra::SiconosMatrix>(1, 1);
  _M->setIdentity();
}

void user_defined::MyDS::computeF(double t) { _f->setValue(0, 0); }
void user_defined::MyDS::computeF(double, std::shared_ptr<siconos::algebra::SiconosVector>) {
  _f->setValue(0, 0);
}


void user_defined::MyDS::computeJacobianfOver_x(
    const Eigen::Ref<siconos::algebra::SiconosVector> &state, double time) {
  jacobianfVectorOver_x_->setValue(0, 0, 0);
}

void user_defined::MyDS::computeRhs(double t) { ; }
void user_defined::MyDS::resetNonSmoothPart(unsigned int level) { _r->setZero(); }
