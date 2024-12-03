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

// #define SICONOS_DEBUG

user_defined::MyDS::MyDS(std::shared_ptr<siconos::algebra::SiconosVector> x0)
    : FirstOrderNonLinearDS(x0) {
  jacobianfVectorOver_x_ = std::make_shared<siconos::algebra::SiconosMatrix>(4, 4);
  _f = std::make_shared<siconos::algebra::SiconosVector>(4);
  _M = std::make_shared<siconos::algebra::SiconosMatrix>(4, 4);
  _M->setIdentity();

  Q = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  Q->setIdentity();
  K1 = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  K1->setValue(0, 0, 0.0);
  K1->setValue(0, 1, 1.0 / 2.0);
  K1->setValue(1, 0, -1.0 / 2.0);
  K1->setValue(1, 1, +1.0);
}

void user_defined::MyDS::computefVector(
    const Eigen::Ref<siconos::algebra::SiconosVector> &state, double time) {
  auto QX = std::make_shared<siconos::algebra::SiconosVector>(2);
  auto X = std::make_shared<siconos::algebra::SiconosVector>(2);

  X->setValue(0, (state->getValue(0) - 2.0));
  X->setValue(1, (state->getValue(1) + 1.0));

  *QX = *Q * *X;

  auto K1P = std::make_shared<siconos::algebra::SiconosVector>(2);
  auto P = std::make_shared<siconos::algebra::SiconosVector>(2);
  P->setValue(0, state->getValue(2));
  P->setValue(1, state->getValue(3));
  *K1P = *K1 * *P;

  auto alphatmp = std::make_shared<siconos::algebra::SiconosVector>(2);

  alpha(time, state, alphatmp);

  _f->setValue(0, alphatmp->getValue(0));
  _f->setValue(1, alphatmp->getValue(1));
  _f->setValue(2, -QX->getValue(0) + K1P->getValue(0));
  _f->setValue(3, -QX->getValue(1) + K1P->getValue(1));
}

void user_defined::MyDS::computeJacobianfOver_x(
    const Eigen::Ref<siconos::algebra::SiconosVector> &state, double time) {
  auto jacXalpha = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);

  JacobianXalpha(t, state, jacXalpha);

  jacobianfVectorOver_x_->setValue(0, 0, jacXalpha->getValue(0, 0));
  jacobianfVectorOver_x_->setValue(0, 1, jacXalpha->getValue(0, 1));
  jacobianfVectorOver_x_->setValue(0, 2, 0.0);
  jacobianfVectorOver_x_->setValue(0, 3, 0.0);
  jacobianfVectorOver_x_->setValue(1, 0, jacXalpha->getValue(1, 0));
  jacobianfVectorOver_x_->setValue(1, 1, jacXalpha->getValue(1, 1));
  jacobianfVectorOver_x_->setValue(1, 2, 0.0);
  jacobianfVectorOver_x_->setValue(1, 3, 0.0);
  jacobianfVectorOver_x_->setValue(2, 0, -Q->getValue(0, 0));
  jacobianfVectorOver_x_->setValue(2, 1, -Q->getValue(0, 1));
  jacobianfVectorOver_x_->setValue(2, 2, K1->getValue(0, 0));
  jacobianfVectorOver_x_->setValue(2, 3, K1->getValue(0, 1));
  jacobianfVectorOver_x_->setValue(3, 0, -Q->getValue(1, 0));
  jacobianfVectorOver_x_->setValue(3, 1, -Q->getValue(1, 1));
  jacobianfVectorOver_x_->setValue(3, 2, K1->getValue(1, 0));
  jacobianfVectorOver_x_->setValue(3, 3, K1->getValue(1, 1));
}

void user_defined::MyDS::alpha(double t,
                               std::shared_ptr<siconos::algebra::SiconosVector> state,
                               std::shared_ptr<siconos::algebra::SiconosVector> _alpha) {
  _alpha->setValue(0, 1.0 / 2.0 * state->getValue(1) + 1.0 / 2.0);
  _alpha->setValue(1, -1.0 / 2.0 * state->getValue(0) - state->getValue(1));
}

void user_defined::MyDS::JacobianXalpha(
    double t, std::shared_ptr<siconos::algebra::SiconosVector> state,
    std::shared_ptr<siconos::algebra::SiconosMatrix> JacXalpha) {
  JacXalpha->setValue(0, 0, 0.0);
  JacXalpha->setValue(0, 1, 1.0 / 2.0);
  JacXalpha->setValue(1, 0, -1.0 / 2.0);
  JacXalpha->setValue(1, 1, -1.0);

#ifdef SICONOS_DEBUG
  std::cout << "JacXalpha\n" << std::endl;
  ;
  JacXalpha->display();
#endif
}
