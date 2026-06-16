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

#ifndef YOYOH
#define YOYOH

#include <array>
#include <cmath>

namespace user {

inline constexpr double r = 0.005;  // rayons des cylindres du yoyo
inline constexpr double R = 0.03;
inline constexpr double m = 0.1;  // masse du yoyo
inline constexpr double I =
    0.5 * m * R * R;              // moment d'inertie du yoyo par rapport à l'axe de rotation
inline constexpr double L = 0.4;  // longueur totale du fil
inline constexpr double g = 9.8;  // valeur de champ de gravité
inline constexpr double epsilon =
    0.000;  // coefficient de frottement entre le tambour et le fil

inline constexpr double e_eq =
    (I - m * r * r) / (I + m * r * r);  // coefficient de restitution équivalent
inline constexpr double k0 = (1 - e_eq) / (1 + e_eq);
inline const double ks = (1 - e_eq) * (3 + e_eq * e_eq) / std::pow(1 + e_eq, 3);
inline const double w = ks + 0.005;  // coefficient d'accélération
inline const double Cy = std::sqrt(2 * L * (I + m * r * r) / (m * r * r * g));
inline constexpr double c1 = 10;
inline constexpr double c2 = 60;  // coefficients d'amortissement

inline constexpr double A[] = {1.369, -3.665, 1.748, -0.4906, 0.01883};

// amplitudes objectifs  avec les instants corrependants
inline constexpr int G = 5;
inline constexpr std::array<double, G> Som = {L / 2, L / 4, L, L, L / 2};

inline constexpr std::array<double, G> times = {5, 20, 35, 50, 60};

inline const double ue = std::pow(e_eq * w + e_eq - 1 + w, 2) /
                         ((std::pow(1 + e_eq, 2) * w + std::pow(1 - e_eq, 2)) * w);

// Fonction thetaset
inline double thetaset(const std::array<double, G>& Som, size_t i) {
  assert(i < G);
  return (1 - ue) * (Som[i]) / r;
}

// inline double accelerationmain(int n, double A[], double Cy, double time) {
//   double resultat = 0;
//   for (int i = 0; i < n; i++) resultat = resultat + A[i] * sin((i + 1) * (M_PI / Cy) *
//   time); return resultat;
// }

}  // namespace user

#endif
