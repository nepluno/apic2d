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

#include "openglutils.h"

#ifdef __APPLE__
#include <GLUT/glut.h>  // why does Apple have to put glut.h here...
#else
#include <GL/glut.h>  // ...when everyone else puts it here?
#endif

#include "math_defs.h"

void draw_circle2d(const Vector2s& centre, scalar rad, int segs) {
  for (int i = 0; i <= segs; i++) {
    scalar cosine = rad * cos(i * 2 * 3.14159 / (scalar)(segs));
    scalar sine = rad * sin(i * 2 * 3.14159 / (scalar)(segs));
    Vector2s tmp = (Vector2s(cosine, sine) + centre);
    glVertex2fv(tmp.data());
    if (i > 0 && i < segs) {
      glVertex2fv(tmp.data());        
    }
  }
}

void draw_grid2d(const Vector2s& origin, scalar dx, int nx, int ny) {
  scalar width = nx * dx;
  scalar height = ny * dx;

  for (int i = 0; i <= nx; i++) {
    Vector2s a(i * dx, 0);
    Vector2s b(i * dx, height);
    Vector2s oa = origin + a;
    Vector2s ob = origin + b;
    glVertex2fv(oa.data());
    glVertex2fv(ob.data());
  }
  for (int j = 0; j <= ny; ++j) {
    Vector2s a(0, j * dx);
    Vector2s b(width, j * dx);
    Vector2s oa = origin + a;
    Vector2s ob = origin + b;
    glVertex2fv(oa.data());
    glVertex2fv(ob.data());
  }
}

void draw_box2d(const Vector2s& origin, scalar width, scalar height) {
  Vector2s o1 = origin + Vector2s(0, height);
  Vector2s o2 = origin + Vector2s(width, height);
  Vector2s o3 = origin + Vector2s(width, 0);

  glVertex2fv(origin.data());
  glVertex2fv(o1.data());
  glVertex2fv(o1.data());
  glVertex2fv(o2.data());
  glVertex2fv(o2.data());
  glVertex2fv(o3.data());
  glVertex2fv(o3.data());
  glVertex2fv(origin.data());
}

void draw_segmentset2d(const std::vector<Vector2s>& vertices, const std::vector<Vector2i>& edges) {
  for (unsigned int i = 0; i < edges.size(); ++i) {
    glVertex2fv(vertices[edges[i][0]].data());
    glVertex2fv(vertices[edges[i][1]].data());
  }
}

void draw_points2d(const std::vector<Vector2s>& points) {
  for (unsigned int i = 0; i < points.size(); ++i) {
    glVertex2fv(points[i].data());
  }
}

void draw_polygon2d(const std::vector<Vector2s>& vertices) {
  for (unsigned int i = 0; i <= vertices.size(); ++i) {
    glVertex2fv(vertices[i % vertices.size()].data());
    if (i > 0 && i < vertices.size()) {
      glVertex2fv(vertices[i].data());    
    }
  }
}

void draw_polygon2d(const std::vector<Vector2s>& vertices, const std::vector<int>& order) {
  for (unsigned int i = 0; i <= order.size(); ++i) {
    glVertex2fv(vertices[order[i % order.size()]].data());
    if (i > 0 && i < order.size()) {
      glVertex2fv(vertices[order[i]].data());    
    }
  }
}
void draw_segment2d(const Vector2s& start, const Vector2s& end) {
  glVertex2fv(start.data());
  glVertex2fv(end.data());
}

void draw_arrow2d(const Vector2s& start, const Vector2s& end, scalar arrow_head_len) {
  Vector2s direction = end - start;

  Vector2s dir_norm = direction;

  // TODO Possibly automatically scale arrowhead length based on vector
  // magnitude
  if (dir_norm.norm() < 1e-14) return;

  dir_norm.normalize();
  Vector2s perp(dir_norm[1], -dir_norm[0]);

  Vector2s tip_left = end + arrow_head_len / (scalar)sqrt(2.0) * (-dir_norm + perp);
  Vector2s tip_right = end + arrow_head_len / (scalar)sqrt(2.0) * (-dir_norm - perp);

  glVertex2fv(start.data());
  glVertex2fv(end.data());
  glVertex2fv(end.data());
  glVertex2fv(tip_left.data());
  glVertex2fv(end.data());
  glVertex2fv(tip_right.data());
}

void draw_trimesh2d(const std::vector<Vector2s>& vertices, const std::vector<Vector3i>& tris) {
  for (unsigned int i = 0; i < tris.size(); ++i) {
    glVertex2fv(vertices[tris[i][0]].data());
    glVertex2fv(vertices[tris[i][1]].data());
    glVertex2fv(vertices[tris[i][2]].data());
  }
}

void hueToRGB(scalar hue, scalar sat, scalar val, scalar& r, scalar& g, scalar& b) {
  // compute hue (adapted from an older Wikipedia article)
  int Hi = (int)(floor(hue / 60.0f)) % 6;
  scalar f = hue / 60 - Hi;
  scalar p = val * (1 - sat);
  scalar q = val * (1 - f * sat);
  scalar t = val * (1 - (1 - f) * sat);

  switch (Hi) {
    case 0:
      r = val;
      g = t;
      b = p;
      break;
    case 1:
      r = q;
      g = val;
      b = p;
      break;
    case 2:
      r = p;
      g = val;
      b = t;
      break;
    case 3:
      r = p;
      g = q;
      b = val;
      break;
    case 4:
      r = t;
      g = p;
      b = val;
      break;
    case 5:
      r = val;
      g = p;
      b = q;
      break;
  }
}

void draw_grid_data2d(Array2s& data, Vector2s origin, scalar dx, bool color) {
  scalar max_val = -1e+37f;
  scalar min_val = 1e+37f;
  for (int j = 0; j < data.nj; ++j)
    for (int i = 0; i < data.ni; ++i) {
      max_val = std::max(data(i, j), max_val);
      min_val = std::min(data(i, j), min_val);
    }

  for (int j = 0; j < data.nj; ++j) {
    for (int i = 0; i < data.ni; ++i) {
      Vector2s bl = origin + Vector2s(i * dx, j * dx);
      scalar r, g, b;
      if (color) {
        hueToRGB(240 * (data(i, j) - min_val) / (max_val - min_val), 1, 1, r, g, b);
      } else {
        scalar gray = (data(i, j) - min_val) / (max_val - min_val);
        r = g = b = gray;
      }
      // TODO Black body colormap, if I can find it.
      glColor3f(r, g, b);
      draw_box2d(bl, dx, dx);
    }
  }
}
