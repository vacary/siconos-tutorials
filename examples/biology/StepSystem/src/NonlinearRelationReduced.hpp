#ifndef NONLINEARRELATIONREDUCED_H
#define NONLINEARRELATIONREDUCED_H

#include <SiconosKernel.hpp>

namespace user_defined {
class NonlinearRelationReduced : public siconos::modeling::FirstOrderType2R {
 protected:
 public:
  virtual ~NonlinearRelationReduced() noexcept = default;

  /** default function to compute h
   *  \param double : current time
   */
  virtual void computeh(double t, const siconos::algebra::BlockVector& x,
                        const siconos::algebra::SiconosVector& lambda,
                        siconos::algebra::SiconosVector& y) override;

  /** default function to compute g
   *  \param double : current time
   */
  virtual void computeg(double t, const siconos::algebra::SiconosVector& lambda,
                        siconos::algebra::BlockVector& r) override;

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJachx(double t, const siconos::algebra::BlockVector& x,
                            const siconos::algebra::SiconosVector& lambda,
                            siconos::algebra::SimpleMatrix& C) override;
  virtual void computeJachlambda(double t, const siconos::algebra::BlockVector& x,
                                 const siconos::algebra::SiconosVector& lambda,
                                 siconos::algebra::SimpleMatrix& D) override;

  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default
   * value .
   */
  virtual void computeJacglambda(double t, const siconos::algebra::SiconosVector& lambda,
                                 siconos::algebra::SimpleMatrix& B) override;
};
}  // namespace user_defined

#endif
