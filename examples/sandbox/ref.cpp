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

/*!\file BouncingBallTS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include <ReferenceClasses.hpp>
#include <SiconosKernel.hpp>
#include <chrono>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

void example() {
  siconos::algebra::SiconosVector q0{7}, q01{7}, velocity0{6};

  siconos::algebra::SiconosMatrix inertia{3, 3};
  // siconos::algebra::SiconosMatrix mass{6, 6};
  auto mass = 10.;
  q0 << 1, 2, 3, 0., 1, 0, 0;
  q01 << 1, 2, 3, 1, 0, 0, 0.;
  velocity0 << 4., 5, 6, 7, 8, 9;
  inertia(0, 0) = 1;
  inertia(1, 1) = 2;
  inertia(2, 2) = 3;

  auto ds = std::make_shared<siconos::modeling::NewtonEulerDS>(q0, velocity0, mass, inertia);

  // auto mass_func = [](Eigen::Ref<siconos::algebra::MapVectorType> pos, double time,
  //                     Eigen::Ref<siconos::algebra::MapType> result) {
  //   result.setZero();
  //   result(0, 0) = 1;
  //   result(1, 1) = 2.;
  //   result(2, 2) = 3.;
  // };
}

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int ndof = 3;  // 100000000;  // degrees of freedom for the ball
    double m = 1;           // Ball mass
    double g = 9.81;        // Gravity

    // -- Initial positions and velocities --
    //          Vector q0{ndof};
    // auto q0 = std::make_shared<Vector>(ndof);
    // q0->setZero();
    //(*q0)(0) = 1.2;

    Vector q0{ndof};
    q0.setZero();
    q0(0) = 1.;

    Vector v0{ndof};
    v0.setZero();
    v0(0) = 1.;
    Vector z0{ndof};
    z0.setZero();
    z0(0) = 1.;

    z0 = v0 + q0;

    // -- Build
    // auto ball = std::make_shared<siconos::internal::devel_model::ClassA>(q0);

    // Vector fext{ndof};
    // fext.setZero();
    // fext(1) = 112.2;
    // ball->setConstantVectorName2(fext);

    // //ball->vectorName2()->display();
    // ball->computeVectorName2(0.4);
    // //ball->vectorName2()->display();

    // //ball->display();

    std::cout << "Second test ... \n";
    // -- Set external forces (weight) --
    auto myforces = [m, g](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
      int i = 0;
      for (auto& v : result) v = time * i++;

      //  siconos::tools::print("call plugin", result);
    };

    auto ball2 = std::make_shared<siconos::internal::devel_model::ClassA>(q0);

    std::cout << (*ball2->vectorName1())(0) << " " << (*ball2->vectorName3())(0) << "\n";

    ball2->computeVectorName2(0.);
    if (ball2->hasVectorName2()) {
      ball2->vectorName2()->display();
    }

    //    ball2->display();
    ball2->setComputeVectorName2Function(myforces);
    // //ball2->display();
    // //ball2->vectorName2()->display();
    ball2->computeVectorName2(0.);
    // //ball2->vectorName2()->display();
    ball2->computeVectorName2(1.);
    ball2->vectorName2()->display();
    if (ball2->hasVectorName2()) {
      ball2->vectorName2()->display();

      auto res = 3 * ball2->vectorName2_view();
      std::cout << res << "\n";
    }
    ball2->vectorName2()->display();

    /// ----- with vectorNameDirect -----

    ball2->computeVectorNameDirect(0.);
    // if (ball2->hasVectorNameDirect()) {
    //   ball2->vectorNameDirect()->display();
    // }

    //    ball2->display();
    ball2->setComputeVectorNameDirectFunction(myforces);
    // //ball2->display();
    // //ball2->vectorNameDirect()->display();
    ball2->computeVectorNameDirect(0.);
    // //ball2->vectorNameDirect()->display();
    ball2->computeVectorNameDirect(1.);
    // ball2->vectorNameDirect()->display();
    if (ball2->hasVectorNameDirect()) {
      //   //   ball2->vectorNameDirect()->display();

      auto res = 3 * ball2->vectorNameDirect_view();
      std::cout << res << "\n";
    }
    // // ball2->vectorNameDirect()->display();
    ////  --------------------------------

    /// ----- with vectorNameSpan -----
    // -- Set external forces (weight) --
    auto myforces_span = [m, g](double time, std::span<double> result) {
      int i = 0;
      for (auto& v : result) v = time * i++;

      result[1] = 12;
      //  siconos::tools::print("call plugin", result);
    };

    ball2->computeVectorNameSpan(0.);
    if (ball2->hasVectorNameSpan()) {
      ball2->vectorNameSpan()->display();
    }

    //    ball2->display();
    ball2->setComputeVectorNameSpanFunction(myforces_span);
    // //ball2->display();
    // //ball2->vectorNameSpan()->display();
    ball2->computeVectorNameSpan(0.);
    // //ball2->vectorNameSpan()->display();
    ball2->computeVectorNameSpan(1.);
    ball2->vectorNameSpan()->display();
    if (ball2->hasVectorNameSpan()) {
      ball2->vectorNameSpan()->display();

      auto res = 3 * ball2->vectorNameSpan_view();
      std::cout << res << "\n";
    }
    ball2->vectorNameSpan()->display();
    ////  --------------------------------

    siconos::algebra::SiconosMatrix mass{ndof, ndof};
    mass.setZero();
    mass(1, 2) = -m * g;
    mass(0, 1) = 12;

    // mass(3,4) = 12; // ça marche avec ndof = 3, pourquoi ????

    ball2->setConstantMatrixName(mass);

    ball2->matrixName()->display();

    auto mass_func = [m, g](Eigen::Ref<siconos::algebra::MapVectorType> pos, double time,
                            Eigen::Ref<siconos::algebra::MapType> result) {
      int i = 0;
      // for (auto& v : result) v = time * i++;

//      result << 1, 2, 3, 4;

      result(2, 2) = -m * g;

      //  siconos::tools::print("call plugin", result);
    };

    ball2->setComputeMatrixNameFunction(mass_func);
    ball2->matrixName()->display();

    Vector pos{ndof};
    pos.setZero();
    pos(1) = 8;
    ball2->computeMatrixName(1., pos);
    ball2->matrixName()->display();

    example();

    return 0;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
