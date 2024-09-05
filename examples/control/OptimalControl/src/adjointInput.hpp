#ifndef ADJOINTINPUT_H
#define ADJOINTINPUT_H

#include <SiconosKernel.hpp>

namespace user_defined {

class adjointInput : public siconos::modeling::FirstOrderNonLinearR {
 protected:
  std::shared_ptr<siconos::algebra::SiconosMatrix> K2{nullptr};

 public:
  adjointInput() = default;
  virtual ~adjointInput() noexcept = default;

  virtual void initialize(siconos::modeling::Interaction &inter) override;

  /** default function to compute h
   *  \param double : current time
   */
  virtual void computeh(double time, const siconos::algebra::BlockVector &x,
                        const siconos::algebra::SiconosVector &lambda,
                        siconos::algebra::BlockVector &z,
                        siconos::algebra::SiconosVector &y) override;

  /** default function to compute g
   *  \param double time, Interaction& inter : current time
   */
  virtual void computeg(double time, const siconos::algebra::BlockVector &x,
                        const siconos::algebra::SiconosVector &lambda,
                        siconos::algebra::BlockVector &z,
                        siconos::algebra::BlockVector &r) override;

  /** default function to compute jacobianH
   *  \param double time, Interaction& inter : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJachx(double time, const siconos::algebra::BlockVector &x,
                            const siconos::algebra::SiconosVector &lambda,
                            siconos::algebra::BlockVector &z,
                            siconos::algebra::SiconosMatrix &C) override;
  virtual void computeJachlambda(double time, const siconos::algebra::BlockVector &x,
                                 const siconos::algebra::SiconosVector &lambda,
                                 siconos::algebra::BlockVector &z,
                                 siconos::algebra::SiconosMatrix &C) override;

  /** default function to compute jacobianG according to lambda
   *  \param double time, Interaction& inter : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default
   * value .
   */
  virtual void computeJacgx(double time, const siconos::algebra::BlockVector &x,
                            const siconos::algebra::SiconosVector &lambda,
                            siconos::algebra::BlockVector &z,
                            siconos::algebra::SiconosMatrix &K) override;
  virtual void computeJacglambda(double time, const siconos::algebra::BlockVector &x,
                                 const siconos::algebra::SiconosVector &lambda,
                                 siconos::algebra::BlockVector &z,
                                 siconos::algebra::SiconosMatrix &B) override;

  double source(double t);

  void beta(double t, const siconos::algebra::BlockVector &xvalue,
            std::shared_ptr<siconos::algebra::SiconosVector> alpha);

  void JacobianXbeta(double t, const siconos::algebra::BlockVector &xvalue,
                     std::shared_ptr<siconos::algebra::SiconosMatrix> JacbetaX);
};
}  // namespace user_defined

#endif
