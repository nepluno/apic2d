// Copyright 2016 Raymond Yun Fei, Christopher Batty, Robert Bridson
//
// Licensed under the Apache License,
// Version 2.0(the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef APIC2D_BLAS_WRAPPER_H_
#define APIC2D_BLAS_WRAPPER_H_

// Simple placeholder code for BLAS calls - replace with calls to a real BLAS
// library

#include <vector>
#include <Eigen/Core>

namespace robertbridson {

namespace BLAS {

// dot products ==============================================================

inline float dot(const std::vector<float> &x, const std::vector<float> &y) {
  // return cblas_ddot((int)x.size(), &x[0], 1, &y[0], 1);

  return Eigen::Map<const Eigen::VectorXf>(x.data(), x.size()).dot(Eigen::Map<const Eigen::VectorXf>(y.data(), y.size()));
}

// inf-norm (maximum absolute value: index of max returned) ==================

inline int index_abs_max(const std::vector<float> &x) {
  // return cblas_idamax((int)x.size(), &x[0], 1);

  int maxind = 0;
  Eigen::Map<const Eigen::VectorXf>(x.data(), x.size()).maxCoeff(&maxind);
  return maxind;
}

// inf-norm (maximum absolute value) =========================================
// technically not part of BLAS, but useful

inline float abs_max(const std::vector<float> &x) { return std::fabs(x[index_abs_max(x)]); }

// saxpy (y=alpha*x+y) =======================================================

inline void add_scaled(float alpha, const std::vector<float> &x, std::vector<float> &y) {
  // cblas_daxpy((int)x.size(), alpha, &x[0], 1, &y[0], 1);
  Eigen::Map<Eigen::VectorXf>(y.data(), y.size()) += Eigen::Map<const Eigen::VectorXf>(x.data(), x.size()) * alpha;
}
}  // namespace BLAS
}  // namespace robertbridson
#endif  // APIC2D_BLAS_WRAPPER_H_
