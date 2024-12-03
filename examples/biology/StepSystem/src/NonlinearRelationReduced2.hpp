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

#ifndef NONLINEARRELATIONREDUCED2_H
#define NONLINEARRELATIONREDUCED2_H

#include <FirstOrderType2R.hpp>

namespace user_defined {
class NonlinearRelationReduced : public siconos::modeling::FirstOrderType2R {
 protected:
 public:
  NonlinearRelationReduced();

  virtual ~NonlinearRelationReduced() noexcept = default;
  void initialize(Interaction& inter) override;
};
}  // namespace user_defined

#endif
