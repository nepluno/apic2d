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

#include <cfloat>

#include "math_defs.h"

void draw_circle2d(const Vector2s& centre, scalar rad, int segs) {
  glBegin(GL_POLYGON);
  for (int i = 0; i < segs; i++) {
    scalar cosine = rad * cos(i * 2 * M_PI / (scalar)(segs));
    scalar sine = rad * sin(i * 2 * M_PI / (scalar)(segs));
    Vector2s tmp = (Vector2s(cosine, sine) + centre);
    glVertex2dv(tmp.data());
  }
  glEnd();
}

void draw_grid2d(const Vector2s& origin, scalar dx, int nx, int ny) {
  scalar width = nx * dx;
  scalar height = ny * dx;

  glBegin(GL_LINES);
  for (int i = 0; i <= nx; i++) {
    Vector2s a(i * dx, 0);
    Vector2s b(i * dx, height);
    Vector2s oa = origin + a;
    Vector2s ob = origin + b;
    glVertex2dv(oa.data());
    glVertex2dv(ob.data());
  }
  for (int j = 0; j <= ny; ++j) {
    Vector2s a(0, j * dx);
    Vector2s b(width, j * dx);
    Vector2s oa = origin + a;
    Vector2s ob = origin + b;
    glVertex2dv(oa.data());
    glVertex2dv(ob.data());
  }
  glEnd();
}

void draw_box2d(const Vector2s& origin, scalar width, scalar height) {
  Vector2s o1 = origin + Vector2s(0, height);
  Vector2s o2 = origin + Vector2s(width, height);
  Vector2s o3 = origin + Vector2s(width, 0);

  glBegin(GL_POLYGON);
  glVertex2dv(origin.data());
  glVertex2dv(o1.data());
  glVertex2dv(o2.data());
  glVertex2dv(o3.data());
  glEnd();
}

void draw_segmentset2d(const std::vector<Vector2s>& vertices,
                       const std::vector<Vector2i>& edges) {
  glBegin(GL_LINES);
  for (unsigned int i = 0; i < edges.size(); ++i) {
    glVertex2dv(vertices[edges[i][0]].data());
    glVertex2dv(vertices[edges[i][1]].data());
  }
  glEnd();
}

void draw_points2d(const std::vector<Vector2s>& points) {
  glBegin(GL_POINTS);
  for (unsigned int i = 0; i < points.size(); ++i) {
    glVertex2dv(points[i].data());
  }
  glEnd();
}

void draw_polygon2d(const std::vector<Vector2s>& vertices) {
  glBegin(GL_POLYGON);
  for (unsigned int i = 0; i < vertices.size(); ++i)
    glVertex2dv(vertices[i].data());
  glEnd();
}

void draw_polygon2d(const std::vector<Vector2s>& vertices,
                    const std::vector<int>& order) {
  glBegin(GL_POLYGON);
  for (unsigned int i = 0; i < order.size(); ++i)
    glVertex2dv(vertices[order[i]].data());
  glEnd();
}
void draw_segment2d(const Vector2s& start, const Vector2s& end) {
  glBegin(GL_LINES);
  glVertex2dv(start.data());
  glVertex2dv(end.data());
  glEnd();
}

void draw_arrow2d(const Vector2s& start, const Vector2s& end,
                  scalar arrow_head_len) {
  Vector2s direction = end - start;

  Vector2s dir_norm = direction;

  // TODO Possibly automatically scale arrowhead length based on vector
  // magnitude
  if (dir_norm.norm() < 1e-14) return;

  dir_norm.normalize();
  Vector2s perp(dir_norm[1], -dir_norm[0]);

  Vector2s tip_left =
      end + arrow_head_len / (scalar)sqrt(2.0) * (-dir_norm + perp);
  Vector2s tip_right =
      end + arrow_head_len / (scalar)sqrt(2.0) * (-dir_norm - perp);

  glBegin(GL_LINES);
  glVertex2dv(start.data());
  glVertex2dv(end.data());
  glVertex2dv(end.data());
  glVertex2dv(tip_left.data());
  glVertex2dv(end.data());
  glVertex2dv(tip_right.data());
  glEnd();
}

void draw_trimesh2d(const std::vector<Vector2s>& vertices,
                    const std::vector<Vector3i>& tris) {
  glBegin(GL_TRIANGLES);
  for (unsigned int i = 0; i < tris.size(); ++i) {
    glVertex2dv(vertices[tris[i][0]].data());
    glVertex2dv(vertices[tris[i][1]].data());
    glVertex2dv(vertices[tris[i][2]].data());
  }
  glEnd();
}

void hueToRGB(scalar hue, scalar sat, scalar val, scalar& r, scalar& g,
              scalar& b) {
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
  scalar max_val = FLT_MIN;
  scalar min_val = FLT_MAX;
  for (int j = 0; j < data.nj; ++j)
    for (int i = 0; i < data.ni; ++i) {
      max_val = std::max(data(i, j), max_val);
      min_val = std::min(data(i, j), min_val);
    }

  for (int j = 0; j < data.nj; ++j) {
    for (int i = 0; i < data.ni; ++i) {
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      Vector2s bl = origin + Vector2s(i * dx, j * dx);
      scalar r, g, b;
      if (color) {
        hueToRGB(240 * (data(i, j) - min_val) / (max_val - min_val), 1, 1, r, g,
                 b);
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

void draw_trimesh3d(const std::vector<Vector3s>& vertices,
                    const std::vector<Vector3i>& tris) {
  glBegin(GL_TRIANGLES);
  for (unsigned int i = 0; i < tris.size(); ++i) {
    glVertex3dv(vertices[tris[i][0]].data());
    glVertex3dv(vertices[tris[i][1]].data());
    glVertex3dv(vertices[tris[i][2]].data());
  }
  glEnd();
}

void draw_trimesh3d(const std::vector<Vector3s>& vertices,
                    const std::vector<Vector3i>& tris,
                    const std::vector<Vector3s>& normals) {
  glBegin(GL_TRIANGLES);
  for (unsigned int i = 0; i < tris.size(); ++i) {
    glNormal3dv(normals[tris[i][0]].data());
    glVertex3dv(vertices[tris[i][0]].data());
    glNormal3dv(normals[tris[i][1]].data());
    glVertex3dv(vertices[tris[i][1]].data());
    glNormal3dv(normals[tris[i][2]].data());
    glVertex3dv(vertices[tris[i][2]].data());
  }
  glEnd();
}

void draw_box3d(const Vector3s& dimensions) {
  // Draw an axis-aligned box with specified dimensions,
  // where the midpoint of the box is at the origin

  scalar width = dimensions[0];
  scalar height = dimensions[1];
  scalar depth = dimensions[2];

  glBegin(GL_POLYGON);
  glNormal3d(-1, 0, 0);
  glVertex3d(-0.5 * width, -0.5 * height, 0.5 * depth);
  glVertex3d(-0.5 * width, 0.5 * height, 0.5 * depth);
  glVertex3d(-0.5 * width, 0.5 * height, -0.5 * depth);
  glVertex3d(-0.5 * width, -0.5 * height, -0.5 * depth);
  glEnd();

  glBegin(GL_POLYGON);
  glNormal3d(1, 0, 0);
  glVertex3d(0.5 * width, -0.5 * height, 0.5 * depth);
  glVertex3d(0.5 * width, 0.5 * height, 0.5 * depth);
  glVertex3d(0.5 * width, 0.5 * height, -0.5 * depth);
  glVertex3d(0.5 * width, -0.5 * height, -0.5 * depth);
  glEnd();

  glBegin(GL_POLYGON);
  glNormal3d(0, 0, -1);
  glVertex3d(-0.5 * width, -0.5 * height, -0.5 * depth);
  glVertex3d(0.5 * width, -0.5 * height, -0.5 * depth);
  glVertex3d(0.5 * width, 0.5 * height, -0.5 * depth);
  glVertex3d(-0.5 * width, 0.5 * height, -0.5 * depth);
  glEnd();

  glBegin(GL_POLYGON);
  glNormal3d(0, 0, 1);
  glVertex3d(-0.5 * width, -0.5 * height, 0.5 * depth);
  glVertex3d(0.5 * width, -0.5 * height, 0.5 * depth);
  glVertex3d(0.5 * width, 0.5 * height, 0.5 * depth);
  glVertex3d(-0.5 * width, 0.5 * height, 0.5 * depth);
  glEnd();

  glBegin(GL_POLYGON);
  glNormal3d(0, -1, 0);
  glVertex3d(-0.5 * width, -0.5 * height, 0.5 * depth);
  glVertex3d(0.5 * width, -0.5 * height, 0.5 * depth);
  glVertex3d(0.5 * width, -0.5 * height, -0.5 * depth);
  glVertex3d(-0.5 * width, -0.5 * height, -0.5 * depth);
  glEnd();

  glBegin(GL_POLYGON);
  glNormal3d(0, 1, 0);
  glVertex3d(-0.5 * width, 0.5 * height, 0.5 * depth);
  glVertex3d(0.5 * width, 0.5 * height, 0.5 * depth);
  glVertex3d(0.5 * width, 0.5 * height, -0.5 * depth);
  glVertex3d(-0.5 * width, 0.5 * height, -0.5 * depth);
  glEnd();
}
