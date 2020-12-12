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

#ifndef APIC2D_OPENGL_UTILS_H_
#define APIC2D_OPENGL_UTILS_H_

#include <vector>

#include "math_defs.h"
#include "array2.h"

void draw_circle2d(const Vector2s& centre, scalar rad, int segs);
void draw_grid2d(const Vector2s& origin, scalar dx, int nx, int ny);
void draw_box2d(const Vector2s& origin, scalar width, scalar height);
void draw_segmentset2d(const std::vector<Vector2s>& vertices,
                       const std::vector<Vector2i>& edges);
void draw_points2d(const std::vector<Vector2s>& points);
void draw_polygon2d(const std::vector<Vector2s>& vertices);
void draw_polygon2d(const std::vector<Vector2s>& vertices,
                    const std::vector<int>& order);
void draw_segment2d(const Vector2s& start, const Vector2s& end);
void draw_arrow2d(const Vector2s& start, const Vector2s& end,
                  scalar arrow_head_len);
void draw_grid_data2d(Array2s& data, Vector2s origin, scalar dx,
                      bool color = false);
void draw_trimesh2d(const std::vector<Vector2s>& vertices,
                    const std::vector<Vector3i>& tris);

void draw_trimesh3d(const std::vector<Vector3s>& vertices,
                    const std::vector<Vector3i>& tris);
void draw_trimesh3d(const std::vector<Vector3s>& vertices,
                    const std::vector<Vector3i>& tris,
                    const std::vector<Vector3s>& normals);
void draw_box3d(const Vector3s& dimensions);

#endif  // APIC2D_OPENGL_UTILS_H_