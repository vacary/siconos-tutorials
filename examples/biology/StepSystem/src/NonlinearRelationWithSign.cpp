#ifndef NONLINEARRELATIONWITHSIGN_CPP
#define NONLINEARRELATIONWITHSIGN_CPP

#include "NonlinearRelationWithSign.hpp"

// #include "const.h"
#define SICONOS_DEBUG

/*y = h(X)*/
void user_defined::NonlinearRelationWithSign::computeh(
    double t, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SiconosVector& y) {
  // siconos::algebra::SiconosVector& lambda = *inter.lambda(0);

#ifdef SICONOS_DEBUG
  std::cout << "******** user_defined::NonlinearRelationWithSign::computeh computeh at " << t
            << std::endl;
#endif

  y.setValue(0, 4.0 - x(0));
  y.setValue(1, 4.0 - x(1));
  y.setValue(2, 8.0 - x(0));
  y.setValue(3, 8.0 - x(1));
#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  y.display();
#endif
}

/*g=g(lambda)*/
void user_defined::NonlinearRelationWithSign::computeg(
    double t, const siconos::algebra::SiconosVector& lambda,
    siconos::algebra::BlockVector& r) {
#ifdef SICONOS_DEBUG
  std::cout << "*** user_defined::NonlinearRelationWithSign::computeg     computeg at: " << t
            << std::endl;
#endif

  r.setValue(0, 10.0 * (1 - lambda(2)) * (1 + lambda(1)));
  r.setValue(1, 10.0 * (1 + lambda(0)) * (1 - lambda(3)));

#ifdef SICONOS_DEBUG
  std::cout << "user_defined::NonlinearRelationWithSign::computeg with lambda=" << std::endl;
  lambda.display();
  std::cout << std::endl;
  std::cout << "user_defined::NonlinearRelationWithSign::computeg modif g_alpha : \n";
  r.display();
  std::cout << std::endl;
#endif
}

void user_defined::NonlinearRelationWithSign::computeJachx(
    double t, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SiconosMatrix& C) {
  C.setValue(0, 0, -1);
  C.setValue(0, 1, 0);
  C.setValue(1, 0, 0);
  C.setValue(1, 1, -1);
  C.setValue(2, 0, -1);
  C.setValue(2, 1, 0);
  C.setValue(3, 0, 0);
  C.setValue(3, 1, -1);

#ifdef SICONOS_DEBUG
  std::cout << "user_defined::NonlinearRelationWithSign::computeJachx computeJachx "
            << " at "
            << " " << t << ":" << std::endl;
  C.display();
  std::cout << std::endl;
#endif
}

void user_defined::NonlinearRelationWithSign::computeJachlambda(
    double t, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SiconosMatrix& D) {
  D.zero();
}

void user_defined::NonlinearRelationWithSign::computeJacglambda(
    double t, const siconos::algebra::SiconosVector& lambda,
    siconos::algebra::SiconosMatrix& B) {
  //  double *g = &(*Jacglambda)(0,0);
  B.setValue(0, 0, 0);
  B.setValue(1, 0, 10.0 * (1 - lambda(3)));
  B.setValue(0, 1, 10.0 * (1 - lambda(2)));
  B.setValue(1, 1, 0);
  B.setValue(0, 2, -10.0 * (1 + lambda(1)));
  B.setValue(1, 2, 0);
  B.setValue(0, 3, 0);
  B.setValue(1, 3, -10.0 * (1 + lambda(0)));

#ifdef SICONOS_DEBUG
  std::cout << "user_defined::NonlinearRelationWithSign::computeJacgx "
            << " at "
            << " " << t << std::endl;
  B.display();
  std::cout << std::endl;
#endif
}
#endif
