#ifndef BALL_RELATION_HPP
#define BALL_RELATION_HPP

#include <SiconosKernel.hpp>

namespace user_defined {

class BallR : public siconos::modeling::LagrangianScleronomousR {
  double alpha = 0.1;

 public:
  void computeh(const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z,
                siconos::algebra::SiconosVector& y) {
    y.setValue(0, q.getValue(0) + alpha * q.getValue(1));
  }

  void computeJacobianhOver_q(const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z) {
    jacobianhOver_q_->setValue(0, 0, 1.0);
    jacobianhOver_q_->setValue(0, 1, alpha);
  }
};
}  // namespace user_defined
#endif
