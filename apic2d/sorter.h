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

#ifndef APIC2D_SORTER_H_
#define APIC2D_SORTER_H_

#include <vector>

#include "math_defs.h"

struct Particle;
class FluidSim;

class sorter {
 public:
  sorter(int ni_, int nj_);
  ~sorter();

  void sort(FluidSim* sim);

  template<typename Callable>
  void getNeigboringParticles_cell(int i, int j, int wl, int wh, int hl, int hh, Callable func) {
    for (int si = i + wl; si <= i + wh; si++)
      for (int sj = j + hl; sj <= j + hh; sj++) {
        if (si < 0 || si > ni - 1 || sj < 0 || sj > nj - 1) continue;
        func(cells[sj * ni + si]);
      }
  }

  int getNumParticleAt(int i, int j);
  void deleteAllParticles();

  std::vector<std::vector<Particle*> > cells;
  int ni;
  int nj;
};

#endif  // APIC2D_SORTER_H_