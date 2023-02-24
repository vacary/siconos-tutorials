

// #ifdef _WIN32
// #define SICONOS_EXPORT extern "C" __declspec(dllexport)
// #else
// #define SICONOS_EXPORT extern "C"
// #endif
#include <stdio.h>

#include <cassert>
#include <cmath>
#include <numbers>

#include "donnees.h"

#undef restrict
#define restrict __restrict

#if defined(__cplusplus)
extern "C" {
#endif

static double puissance(double a, int n) {
  double resultat = 1;

  for (int i = 0; i < n; i++) resultat = resultat * a;
  return resultat;
}

static double ks = (1 - e_eq) * (3 + e_eq * e_eq) / puissance(1 + e_eq, 3);
static double w = ks + 0.005;  // coefficient d'accélération
static double ue =
    pow(e_eq * w + e_eq - 1 + w, 2) / ((pow(1 + e_eq, 2) * w + pow(1 - e_eq, 2)) * w);
static double Cy = sqrt(2 * L * (I + m * r * r) / (m * r * r * g));

static double accelerationmain(int n, double A[], double Cy, double time) {
  double resultat = 0;
  for (int i = 0; i < n; i++) resultat = resultat + A[i] * sin((i + 1) * (M_PI / Cy) * time);
  return resultat;
}

static double vitessemain(int n, double A[], double Cy, double time) {
  double resultat = 0;
  for (int i = 0; i < n; i++)
    resultat = resultat - A[i] * cos((i + 1) * (M_PI / Cy) * time) / ((i + 1) * M_PI / Cy);
  return resultat;
}

// // amplitudes objectifs  avec les instants corrependants
// int const G = 5;
// double temps[G] = {5, 20, 35, 50, 60};
// double Som[G] = {L / 2, L / 4, L, L, L / 2};

static double thetaset(double Som[], size_t i) {
  assert(i < G);
  return (1 - ue) * (Som[i]) / r;
}

// forces exterieures appliquees sur le yoyo dans la phase contrainte
void force_ext(double time, unsigned int sizeOfq, double* restrict fExt, unsigned int sizeZ,
               double* restrict z) {
  fExt[0] = -m * r * g;
  fExt[1] = 0;
  // fExt[2] = accelerationmain(5,A,Cy,time);
  fExt[2] = 0;
}

// forces exterieures appliquees sur le yoyo dans la phase libre
void force_extf(double time, unsigned int sizeOfq, double* restrict fExt, unsigned int sizeZ,
                double* restrict z) {
  fExt[0] = 0;
  fExt[1] = -m * g;
  // fExt[2] = accelerationmain (5,A,Cy,time);
  fExt[2] = 0;
}

// forces interieures appliquees sur le yoyo dans la phase contrainte
void F_int(double time, unsigned int sizeOfq, double* restrict q, double* restrict velocity,
           double* restrict fInt, unsigned int sizeZ, double* restrict z) {
  fInt[0] = r * epsilon * (velocity[0]);
  fInt[1] = 0;
  int i = 0;
  while (temps[i] < time) i++;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i))
    fInt[2] = -w * g;
  else
    fInt[2] = c1 * velocity[2] + c2 * q[2];
  // fInt[2] =0;
}

void jacobianFIntq(double time, unsigned int sizeOfq, double* restrict q,
                   double* restrict velocity, double* restrict jacob, unsigned int sizeZ,
                   double* restrict z) {
  jacob[0] = 0;
  jacob[1] = 0;
  int i = 0;
  while (temps[i] < time) i++;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i))
    jacob[2] = 0;
  else
    jacob[2] = c2;
  // jacob[2] =0;
}

void jacobianVFInt(double time, unsigned int sizeOfq, double* restrict q,
                   double* restrict velocity, double* restrict jacob, unsigned int sizeZ,
                   double* restrict z) {
  jacob[0] = r * epsilon;
  jacob[1] = 0;
  int i = 0;
  while (temps[i] < time) i++;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i))
    jacob[2] = 0;
  else
    jacob[2] = c1;
  //  jacob[2] =0;
}

void h1(unsigned int sizeDS, double* restrict q, double time, unsigned int sizeY,
        double* restrict y, unsigned int sizeOfZ, double* restrict z) {
  y[0] = q[1] - r * q[0] + L - q[2];
}

void G10(unsigned int sizeDS, double* restrict q, double time, unsigned int sizeY,
         double* restrict G, unsigned int sizeOfZ, double* restrict z) {
  G[0] = -r;
  G[1] = 1;
  G[2] = -1;
}

void G11(unsigned int sizeDS, double* restrict q, double time, unsigned int sizeY,
         double* restrict G, unsigned int sizeOfZ, double* restrict z) {
  G[0] = 0;
  // G[0]= -vitessemain(5,A,Cy,time);
}

// forces interieures appliquees sur le yoyo dans la phase libre
void F_intf(double time, unsigned int sizeOfq, double* restrict q, double* restrict velocity,
            double* restrict fInt, unsigned int sizeZ, double* restrict z) {
  fInt[0] = 0;
  fInt[1] = 0;
  int i = 0;
  while (temps[i] < time) i++;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i))
    fInt[2] = -w * g;
  else
    fInt[2] = c1 * velocity[2] + c2 * q[2];
  // fInt[2] =0;
}

void jacobianFIntqf(double time, unsigned int sizeOfq, double* restrict q,
                    double* restrict velocity, double* restrict jacob, unsigned int sizeZ,
                    double* restrict z) {
  jacob[0] = 0;
  jacob[1] = 0;
  int i = 0;
  while (temps[i] < time) i++;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i))
    jacob[2] = 0;
  else
    jacob[2] = c2;
  // jacob[2] =0;
}

void jacobianVFIntf(double time, unsigned int sizeOfq, double* restrict q,
                    double* restrict velocity, double* restrict jacob, unsigned int sizeZ,
                    double* restrict z) {
  jacob[0] = 0;
  jacob[1] = 0;
  int i = 0;
  while (temps[i] < time) i++;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i))
    jacob[2] = 0;
  else
    jacob[2] = c1;
  // jacob[2] =0;
}

#if defined(__cplusplus)
}
#endif
