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

#ifndef MYDSDS_H
#define MYDSDS_H

#include <SiconosKernel.hpp>

namespace user_defined {

class MyDS : public siconos::modeling::FirstOrderNonLinearDS {
 protected:
  std::shared_ptr<siconos::algebra::SiconosMatrix> Q{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> K1{nullptr};
  // SiconosMatrix * K1T;

 public:
  /** default constructor
   * \param initial conditions
   */
  MyDS(Eigen::Ref<siconos::algebra::SiconosVector> x0);

  // ===== DESTRUCTOR =====

  /** destructor
   */
  virtual ~MyDS() noexcept = default;
};

siconos::algebra::SiconosVector alpha(
    double t, const Eigen::Ref<const siconos::algebra::SiconosVector>& xvalue);

siconos::algebra::SiconosMatrix JacobianXalpha(
    double t, const Eigen::Ref<const siconos::algebra::SiconosVector>& xvalue);
}  // namespace user_defined

#endif
