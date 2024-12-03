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

#include <SiconosKernel.hpp>
#include <cassert>
#include <chrono>
#include <cmath>
#include <numbers>

#include "donnees.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

int main(int argc, char* argv[]) {
  try {
    unsigned int nDof = 3;  // nombre de degrés de liberté du yoyo
    double t0 = 0;          //  instants initial et final de la simulation
    double T = 50;
    double h = 0.001;          // pas de discrétisation du temps
    const double theta = 0.5;  // coefficient pour le générateur de simulation
    double e = 0;              // coefficient de restitution
    int N = ceil((T - t0) / h) + 1;

    // déclarations

    auto M = std::make_shared<Matrix>(nDof, nDof);
    Vector q0{nDof};
    q0.setZero();
    Vector v0{nDof};
    v0.setZero();

    auto loi = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianRheonomousR>(
        "YoyoPlugin:h1", "YoyoPlugin:G10", "YoyoPlugin:G11");

    auto H = std::make_shared<Matrix>(1, nDof);
    (*H)(0, 0) = 0;
    (*H)(0, 1) = 0;
    (*H)(0, 2) = 0;
    auto loi0 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation0 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H);

    // pramètres du solveur siconos

    unsigned int outputSize = 9;
    Matrix dataPlot(N, outputSize);

   q0(0) = L / (2 * r);        //  vlaeur de  teta( 0 )
   q0(1) = -L + r *q0(0);  // y ( 0 )
   q0(2) = 0;                  // valeur initiale de h
    // vitesses initiels de teta , y et  h
    v0(0) = 0;
    v0(1) = r * v0(0);
    v0(2) = 0;

    // Objectifs du joueurs avec les instants correspendants
    Matrix Controle(G + 1, 2);
    Controle(0, 1) =q0(1) -q0(2);
    for (int i = 0; i < G; i++) {
      Controle(i + 1, 0) = temps[i];
      Controle(i + 1, 1) = Som[i] - L;
    }

    int k = 0;

    std::shared_ptr<Vector> q, v;
    auto start = std::chrono::system_clock::now();
    while (k < N) {
      ///////////////////////////////////////Phase contrainte
      /////////////////////////////////////////

      M->setIdentity();
      (*M)(0, 0) = I + m * r * r;
      (*M)(0, 2) = m * r;
      (*M)(1, 0) = -r;
      (*M)(1, 1) = 1;
      (*M)(1, 2) = -1;

      if (k != 0) {
        t0 = dataPlot(k - 1, 0);
       q0(0) = (*q)(0);
       q0(2) = (*q)(2);
       q0(1) =q0(2) - L + r *q0(0);
        v0(0) = (*v)(0);
        v0(2) = (*v)(2);
        v0(1) = r * v0(0) + v0(2);
      }

      // création et insertion  du système dynamique représentant la yoyo dans le récipient
      // allDS
      auto yoyo = std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, M);

      yoyo->setComputeFextFunction("YoyoPlugin", "force_ext");
      yoyo->setComputeFIntFunction("YoyoPlugin", "F_int");
      yoyo->setComputeJacobianFIntqDotFunction("YoyoPlugin", "jacobianVFInt");
      yoyo->setComputeJacobianFIntqFunction("YoyoPlugin", "jacobianFIntq");

      ////////////////  loi d'impact et relations /////////////////////////////////

      auto inter = std::make_shared<siconos::modeling::Interaction>(loi0, relation0);

      /////////////////////////  MODEL //////////////////////////////////////////////////
      auto jeu = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
      jeu->insertDynamicalSystem(yoyo);
      jeu->link(inter, yoyo);
      ///////////////////// SIMULATION /////////////////////////////////

      // déscrétisation du temps
      auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

      // -- OneStepIntegrators --
      auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

      // -- OneStepNsProblem --
      auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

      auto s = std::make_shared<siconos::simulation::TimeStepping>(jeu, t, OSI, osnspb);

      q = yoyo->q();
      v = yoyo->velocity();
      auto p = yoyo->p(1);
      auto lambda = inter->lambda(1);

      // --- sauver les valeurs dans une matrice dataPlot

      dataPlot(k, 0) = jeu->t0();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q)(1);
      dataPlot(k, 6) = (*v)(1);
      dataPlot(k, 7) = L + (*q)(1) - r * (*q)(0) - (*q)(2);
      dataPlot(k, 8) = (*q)(1) - (*q)(2);
      k++;

      while (s->hasNextEvent() && (*q)(0) > 0.0) {
        s->computeOneStep();
        // --- Get values to be plotted ---
        dataPlot(k, 0) = s->nextTime();
        dataPlot(k, 1) = (*q)(0);
        dataPlot(k, 2) = (*v)(0);
        dataPlot(k, 3) = (*p)(0);
        dataPlot(k, 4) = (*lambda)(0);  // multiplicateur de lagrange
        dataPlot(k, 5) = (*q)(1);
        dataPlot(k, 6) = (*v)(1);
        dataPlot(k, 7) = L + (*q)(1) - r * (*q)(0) - (*q)(2);  // contrainte géométrique
        dataPlot(k, 8) = (*q)(1) - (*q)(2);
        if (abs((*v)(0)) <= 0.05) cout << "valeur max de theta est : " << (*q)(0) << endl;
        k++;
        s->nextStep();
      }

      //////////////////////////////////////////Phase
      /// libre//////////////////////////////////////////////////////

      M->setIdentity();
      (*M)(1, 1) = m;
      (*M)(0, 0) = I;

      t0 = dataPlot(k - 1, 0) + h;
      if (t0 + h < T) {
       q0(0) = 0;
       q0(2) = (*q)(2);
       q0(1) =q0(2) - L;
        v0(0) = -(*v)(0);
        v0(2) = (*v)(2);
        v0(1) = -r * v0(0) + v0(2);

        yoyo = std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, M);
        yoyo->setComputeFextFunction("YoyoPlugin", "force_extf");
        yoyo->setComputeFIntFunction("YoyoPlugin", "F_intf");
        yoyo->setComputeJacobianFIntqDotFunction("YoyoPlugin", "jacobianVFIntf");
        yoyo->setComputeJacobianFIntqFunction("YoyoPlugin", "jacobianFIntqf");

        inter = std::make_shared<siconos::modeling::Interaction>(loi, relation);

        jeu = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
        jeu->insertDynamicalSystem(yoyo);
        jeu->link(inter, yoyo);

        t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
        OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
        osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();
        s = std::make_shared<siconos::simulation::TimeStepping>(jeu, t, OSI, osnspb);

        q = yoyo->q();
        v = yoyo->velocity();
        p = yoyo->p(1);
        lambda = inter->lambda(1);

        while (s->hasNextEvent()) {
          s->computeOneStep();
          dataPlot(k, 0) = s->nextTime();
          dataPlot(k, 1) = (*q)(0);
          dataPlot(k, 2) = (*v)(0);
          dataPlot(k, 3) = (*p)(0);
          dataPlot(k, 4) = (*lambda)(0);
          dataPlot(k, 5) = (*q)(1);
          dataPlot(k, 6) = (*v)(1);
          dataPlot(k, 7) = L + (*q)(1) - r * (*q)(0) - (*q)(2);
          dataPlot(k, 8) = (*q)(1) - (*q)(2);
          k++;
          if ((*lambda)(0) > 0 && (-r * (*v)(0) + (*v)(1) - (*v)(2)) < 10e-14) break;
          s->nextStep();
        }
      }
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    siconos::algebra::io::write("yoyo.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("fichier.dat", Controle, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "yoyo.ref", eps)) > eps)
      return 1;

    return 0;
    // --- Libérer de la mémoire
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
