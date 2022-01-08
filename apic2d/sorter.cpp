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
  const int np = (int)sim->particles_.size();
  for (int n = 0; n < np; n++) {
    Particle* p = &sim->particles_[n];

    int pi = (int)((p->x_(0) - sim->origin_(0)) / sim->dx_);
    int pj = (int)((p->x_(1) - sim->origin_(1)) / sim->dx_);
    int i = max(0, min(ni - 1, pi));
    int j = max(0, min(nj - 1, pj));
    cells[j * ni + i].push_back(p);
  }
}

int sorter::getNumParticleAt(int i, int j) { return (int)cells[j * ni + i].size(); }

void sorter::deleteAllParticles() {
  for (int j = 0; j < nj; ++j)
    for (int i = 0; i < ni; ++i) {
      cells[j * ni + i].clear();
    }
}