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

static double r = 0.005;  // rayons des cylindres du yoyo
static double R = 0.03;
static double m = 0.1;  // masse du yoyo
static double I =
    0.5 * m * R * R;               // moment d'inertie du yoyo par rapport à l'axe de rotation
static double L = 0.4;          // longueur totale du fil
static double g = 9.8;          // valeur de champ de gravité
static double epsilon = 0.000;  // coefficient de frottement entre le tambour et le fil

static double e_eq = (I - m * r * r) / (I + m * r * r);  // coefficient de restitution équivalent
static double k0 = (1 - e_eq) / (1 + e_eq);
static double c1 = 10;
static double c2 = 60;  // coefficients d'amortissement

static double A[] = {1.369, -3.665, 1.748, -0.4906, 0.01883};

// amplitudes objectifs  avec les instants corrependants
const int G = 5;
static double temps[G] = {5, 20, 35, 50, 60};
static double Som[G] = {L / 2, L / 4, L, L, L / 2};
#endif
