#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <SiconosKernel.hpp>
#include <SiconosMatrix.hpp>
#include <SiconosVector.hpp>
#include <iostream>

namespace py = pybind11;

void computeMassDense(const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
                      Eigen::Ref<siconos::algebra::SiconosDenseMatrix> mass) {
  // mass.setZero();
  mass(0, 0) = 1.;
  mass(1, 1) = 1.;
  mass(2, 2) = 2. / 5. * 0.1 * 0.1;
  // std::cout << "computeM from addon final\n";
}

void computeMassSparse(const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
                       siconos::algebra::SiconosSparseMatrix& mass) {
  mass.setZero();
  for (int i = 0; i < q.size(); ++i) mass.insert(i, i) = 1.0 + q[i] * 0.1;
  mass.makeCompressed();
}

PYBIND11_MODULE(computeM, addons) {
  addons.def("computeMassDense", &computeMassDense);
  addons.def("computeMassSparse", &computeMassSparse);
}
