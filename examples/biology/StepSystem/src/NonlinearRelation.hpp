#ifndef NONLINEARRELATION_H
#define NONLINEARRELATION_H

#include <SiconosKernel.hpp>

namespace user_defined {
class NonlinearRelation : public siconos::modeling::FirstOrderType2R {
 protected:
 public:
  virtual ~NonlinearRelation() noexcept = default;

  /** default function to compute h */
  virtual void computeh(double t, const siconos::algebra::BlockVector& x,
                        const siconos::algebra::SiconosVector& lambda,
                        siconos::algebra::SiconosVector& y) override;

  /** default function to compute g */
  virtual void computeg(double t, const siconos::algebra::SiconosVector& lambda,
                        siconos::algebra::BlockVector& r) override;

  /** default function to compute jacobian of h w.r.t x and lambda */
  virtual void computeJachx(double t, const siconos::algebra::BlockVector& x,
                            const siconos::algebra::SiconosVector& lambda,
                            siconos::algebra::SimpleMatrix& C) override;
  virtual void computeJachlambda(double t, const siconos::algebra::BlockVector& x,
                                 const siconos::algebra::SiconosVector& lambda,
                                 siconos::algebra::SimpleMatrix& D) override;

  /** default function to compute jacobian of g  w.r.t  lambda  */
  virtual void computeJacglambda(double t, const siconos::algebra::SiconosVector& lambda,
                                 siconos::algebra::SimpleMatrix& B) override;
};
}  // namespace user_defined

#endif
