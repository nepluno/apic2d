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

#ifndef APIC2D_KERNEL_H_
#define APIC2D_KERNEL_H_

#include "math_defs.h"

namespace kernel {
inline scalar smooth_kernel(const scalar& r2, const scalar& h) { return std::max(pow(1.0 - r2 / (h * h), 3.0), 0.0); }

inline scalar smooth_kernel_laplacian(const scalar& r2, const scalar& h) {
  scalar x2 = sqrt(r2 / (h * h));
  return x2 > 1.0 ? 0.0 : (1.0 - x2);
}

inline scalar sharp_kernel(const scalar& r2, const scalar& h) { return std::max(h * h / std::max(r2, 1.0e-5f) - 1.0f, 0.0f); }

inline scalar linear_kernel(const Vector2s& d, const scalar& h) { return std::max((1.0 - fabs(d(0) / h)) * (1.0 - fabs(d(1) / h)), 0.0); }

inline scalar quadratic_kernel_1d(const scalar& d, const scalar& h) {
  scalar r = fabs(d) / h;
  if (r < 0.5) {
    return 0.75 - r * r;
  } else if (r < 1.5) {
    return 0.5 * (1.5 - r) * (1.5 - r);
  } else {
    return 0.0;
  }
}

inline scalar quadratic_kernel(const Vector2s& d, const scalar& h) { return quadratic_kernel_1d(d(0), h) * quadratic_kernel_1d(d(1), h); }
}  // namespace kernel

#endif  // APIC2D_KERNEL_H_