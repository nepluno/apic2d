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

#include <cfloat>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "array2_utils.h"
#include "fluidsim.h"
#include "gluvi.h"
#include "openglutils.h"

using namespace std;

// Try changing the grid resolution
int grid_resolution = 100;
scalar timestep = 0.005;
scalar grid_width = 100.0;

FluidSim sim;

// Gluvi stuff
//-------------
Gluvi::PanZoom2D cam(-0.1, -0.35, 1.2);
double oldmousetime;
Vector2s oldmouse;
void display();
void mouse(int button, int state, int x, int y);
void drag(int x, int y);
void timer(int junk);

// Boundary definition - several circles in a circular domain.

Vector2s c0(50, 50), c1(70, 50), c2(30, 35), c3(50, 70);
Vector2s s0(10, 5);
scalar rad0 = 40, rad1 = 10, rad2 = 10, rad3 = 10;
Vector2s o0(0.0, 0.0);

// Main testing code
//-------------
int main(int argc, char **argv) {
  // Setup viewer stuff
  Gluvi::init("Basic Fluid Solver with Static Variational Boundaries", &argc,
              argv);
  Gluvi::camera = &cam;
  Gluvi::userDisplayFunc = display;
  glClearColor(1, 1, 1, 1);

  glutTimerFunc(1000, timer, 0);

  // Set up the simulation
  sim.initialize(o0, grid_width, grid_resolution, grid_resolution, 1.0);

  sim.root_boundary = new FluidSim::Boundary(c0, Vector2s(rad0, 0.0),
                                             FluidSim::BT_CIRCLE, true);

  sim.root_sources = NULL;

  sim.update_boundary();
  sim.init_random_particles();

  Gluvi::run();

  delete sim.root_boundary;

  return 0;
}

void display(void) { sim.render(); }

void timer(int junk) {
  sim.advance(timestep);

  glutPostRedisplay();
  glutTimerFunc(30, timer, 0);
}
