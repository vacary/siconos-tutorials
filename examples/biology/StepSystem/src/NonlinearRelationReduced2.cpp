#ifndef NONLINEARRELATIONREDUCED2_CPP
#define NONLINEARRELATIONREDUCED2_CPP

#include "NonlinearRelationReduced2.hpp"

// #include "const.h"
#define SICONOS_DEBUG

/*y = h(X)*/
void user_defined::NonlinearRelationReduced2::computeh(
    double t, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SiconosVector& y) {
  y(0) = 8.0 - x(0);
  y(1) = 8.0 - x(1);
}

/*g=g(lambda)*/
void user_defined::NonlinearRelationReduced2::computeg(
    double t, const siconos::algebra::SiconosVector& lambda,
    siconos::algebra::BlockVector& r) {
  r.setValue(0, 40.0 * (1 - lambda(0)));
  r.setValue(1, 40.0 * (1 - lambda(1)));
}

void user_defined::NonlinearRelationReduced2::computeJachlambda(
    double t, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SimpleMatrix& D) {
  D.zero();
}

void user_defined::NonlinearRelationReduced2::computeJachx(
    double t, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SimpleMatrix& C) {
  C.setValue(0, 0, -1);
  C.setValue(0, 1, 0);
  C.setValue(1, 0, 0);
  C.setValue(1, 1, -1);
}

void user_defined::NonlinearRelationReduced2::computeJacglambda(
    double t, const siconos::algebra::SiconosVector& lambda,
    siconos::algebra::SimpleMatrix& B) {
  B.setValue(0, 0, -40.0);
  B.setValue(1, 0, 0.0);
  B.setValue(0, 1, 0.0);
  B.setValue(1, 1, -40.0);
}
#endif
