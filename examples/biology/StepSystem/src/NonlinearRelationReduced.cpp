#ifndef NONLINEARRELATIONREDUCED_CPP
#define NONLINEARRELATIONREDUCED_CPP

#include "NonlinearRelationReduced.hpp"

// #include "const.h"
#define SICONOS_DEBUG

/*y = h(X)*/
void user_defined::NonlinearRelationReduced::computeh(
    double t, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SiconosVector& y) {
  y(0) = 8.0 - x(1);
}

/** default function to compute jacobianH
 *  \param double : current time
 *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
 */

/** default function to compute jacobianG according to lambda
 *  \param double : current time
 *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default
 * value .
 */

/*g=g(lambda)*/
void user_defined::NonlinearRelationReduced::computeg(
    double t, const siconos::algebra::SiconosVector& lambda,
    siconos::algebra::BlockVector& r) {
  r.setValue(1, 40.0 * (1 - lambda(0)));
  r.setValue(0, 0.0);
}

void user_defined::NonlinearRelationReduced::computeJachlambda(
    double t, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SiconosMatrix& D) {
  D.zero();
}

void user_defined::NonlinearRelationReduced::computeJachx(
    double t, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SiconosMatrix& C) {
  C(0, 0) = 0;
  C(0, 1) = -1;
}

void user_defined::NonlinearRelationReduced::computeJacglambda(
    double t, const siconos::algebra::SiconosVector& lambda,
    siconos::algebra::SiconosMatrix& B) {
  B(0, 0) = 0.0;
  B(1, 0) = -40.0;
}
#endif
