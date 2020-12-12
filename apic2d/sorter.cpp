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

#include "sorter.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fluidsim.h"

using namespace std;

sorter::sorter(int ni_, int nj_) : ni(ni_), nj(nj_) { cells.resize(ni_ * nj_); }

sorter::~sorter() {}

void sorter::sort(FluidSim* sim) {
  // Clear All Cells
  for (int j = 0; j < nj; ++j)
    for (int i = 0; i < ni; ++i) {
      cells[j * ni + i].clear();
    }

  // Store Into The Cells
  const int np = (int)sim->particles.size();
  for (int n = 0; n < np; n++) {
    Particle* p = &sim->particles[n];

    int pi = (int)((p->x(0) - sim->origin(0)) / sim->dx);
    int pj = (int)((p->x(1) - sim->origin(1)) / sim->dx);
    int i = max(0, min(ni - 1, pi));
    int j = max(0, min(nj - 1, pj));
    cells[j * ni + i].push_back(p);
  }
}

void sorter::getNeigboringParticles_cell(int i, int j, int wl, int wh, int hl,
                                         int hh, std::vector<Particle*>& res) {
  for (int si = i + wl; si <= i + wh; si++)
    for (int sj = j + hl; sj <= j + hh; sj++) {
      if (si < 0 || si > ni - 1 || sj < 0 || sj > nj - 1) continue;
      res.insert(res.end(), cells[sj * ni + si].begin(),
                 cells[sj * ni + si].end());
    }
}

int sorter::getNumParticleAt(int i, int j) {
  return (int)cells[j * ni + i].size();
}

void sorter::deleteAllParticles() {
  for (int j = 0; j < nj; ++j)
    for (int i = 0; i < ni; ++i) {
      cells[j * ni + i].clear();
    }
}