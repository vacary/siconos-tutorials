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

#include "NonlinearRelationReduced2.hpp"

user_defined::NonlinearRelation::NonlinearRelationReduced2() : FirstOrderType2R{} {
  setComputehFunction([](const siconos::algebra::BlockVector& state,
                         const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
                         Eigen::Ref<siconos::algebra::SiconosVector> y) {
    y(0) = 8.0 - x(0);
    y(1) = 8.0 - x(1);
  });

  setComputegFunction([](const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
                         siconos::algebra::BlockVector& res) {
    res.setValue(0, 40.0 * (1 - lambda(0)));
    res.setValue(1, 40.0 * (1 - lambda(1)));
  });
}

void user_defined::NonlinearRelation::NonlinearRelationReduced2::initialize(
    Interaction& inter) {
  FirstOrderR::initialize(inter);

  auto sizeY = inter.dimension();
  auto sizeX = inter.getSizeOfDS();
  auto& DSlink = inter.linkToDSVariables();
  auto sizeZ = DSlink[FirstOrderR::z]->size();

  hasConstantJacobianhOver_state_ = true;
  computejacobianhOver_state_ = nullptr;

  if (!jacobianhOver_state_internal_storage_) {
    jacobianhOver_state_internal_storage_ =
        std::make_unique<std::vector<double>>(sizeY * sizeX);
  }
  jacobianhOver_state_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianhOver_state_internal_storage_->data(), sizeY, sizeX);

  jacobianhOver_state_view_->setZero();
  jacobianhOver_state_view_.setValue(0, 0, -1);
  jacobianhOver_state_view_.setValue(0, 1, 0);
  jacobianhOver_state_view_.setValue(1, 0, 0);
  jacobianhOver_state_view_.setValue(1, 1, -1);

  if (computejacobianhOver_lambda_) {
    if (!jacobianhOver_lambda_internal_storage_) {
      jacobianhOver_lambda_internal_storage_ =
          std::make_unique<std::vector<double>>(sizeY * sizeY);
    }
    jacobianhOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
        jacobianhOver_lambda_internal_storage_->data(), sizeY, sizeY);
    //   relationMat[FirstOrderR::mat_D] =
    // std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeY);
  }

  if (!jacobiangOver_lambda_internal_storage_) {
    jacobiangOver_lambda_internal_storage_ =
        std::make_unique<std::vector<double>>(sizeX * sizeY);
  }
  jacobiangOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobiangOver_lambda_internal_storage_->data(), sizeX, sizeY);

  hasConstantJacobiangOver_lambda_ = true;
  computejacobiangOver_lambda_ = nullptr;

  jacobiangOver_lambda_view_->setZero();
  jacobiangOver_lambda_view_.setValue(0, 0, -40.0);
  jacobiangOver_lambda_view_.setValue(1, 0, 0.0);
  jacobiangOver_lambda_view_.setValue(0, 1, 0.0);
  jacobiangOver_lambda_view_.setValue(1, 1, -40.0);

  checkSize(inter);
}
