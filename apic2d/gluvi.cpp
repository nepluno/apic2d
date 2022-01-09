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

#include "gluvi.h"

#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <fstream>

#include "math_defs.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

namespace Gluvi {
PanZoom2D::PanZoom2D(scalar bottom_, scalar left_, scalar height_) : bottom(bottom_), left(left_), height(height_), action_mode(INACTIVE) {
  default_bottom = bottom;
  default_left = left;
  default_height = height;
}

void PanZoom2D::click(int button, int state, int x, int y) {
  if (state == GLUT_UP) {
    scalar r = height / winheight;
    switch (action_mode) {
      case PAN:
        if (!moved_since_mouse_down) {
          // make mouse click the centre of the window
          left += r * (x - 0.5 * winwidth);
          bottom += r * (0.5 * winheight - y);
          glutPostRedisplay();
        }
        break;
      case ZOOM_IN:
        if (moved_since_mouse_down) {
          // zoom in to selection
          scalar desired_width = fabs((x - clickx) * height / winheight);
          scalar desired_height = fabs((y - clicky) * height / winheight);
          if (desired_height == 0) desired_height = height / winheight;
          if (desired_width * winheight > desired_height * winwidth)
            desired_height = winheight * desired_width / winwidth;
          else
            desired_width = winwidth * desired_height / winheight;
          left += 0.5 * (x + clickx) * height / winheight - 0.5 * desired_width;
          bottom += (winheight - 0.5 * (y + clicky)) * height / winheight - 0.5 * desired_height;
          height = desired_height;
        } else {
          // zoom in by some constant factor on the mouse click
          scalar factor = 0.70710678118654752440084;
          left += (1 - factor) * height * (x / (scalar)winheight);
          bottom += (1 - factor) * height * (1 - y / (scalar)winheight);
          height *= factor;
        }
        glutPostRedisplay();
        break;
      case ZOOM_OUT:
        // zoom out by some constant factor
        {
          scalar factor = 1.41421356237309504880168;
          left -= 0.5 * (factor - 1) * winwidth * height / winheight;
          bottom -= 0.5 * (factor - 1) * height;
          height *= factor;
        }
        glutPostRedisplay();
        break;
      default:;  // nothing to do
    }
    action_mode = INACTIVE;

  } else if (button == GLUT_LEFT_BUTTON)
    action_mode = PAN;
  else if (button == GLUT_MIDDLE_BUTTON) {
    clickx = x;
    clicky = y;
    action_mode = ZOOM_IN;
  } else if (button == GLUT_RIGHT_BUTTON)
    action_mode = ZOOM_OUT;
  moved_since_mouse_down = false;
  oldmousex = x;
  oldmousey = y;
}

void PanZoom2D::drag(int x, int y) {
  if (x != oldmousex || y != oldmousey) {
    moved_since_mouse_down = true;
    if (action_mode == PAN) {
      scalar r = height / winheight;
      left -= r * (x - oldmousex);
      bottom += r * (y - oldmousey);
      glutPostRedisplay();
    } else if (action_mode == ZOOM_IN)
      glutPostRedisplay();
    oldmousex = x;
    oldmousey = y;
  }
}

void PanZoom2D::return_to_default(void) {
  bottom = default_bottom;
  left = default_left;
  height = default_height;
}

void PanZoom2D::transform_mouse(int x, int y, scalar coords[2]) {
  scalar r = height / winheight;
  coords[0] = x * r + left;
  coords[1] = (winheight - y) * r + bottom;
}

void PanZoom2D::gl_transform(void) {
  glViewport(0, 0, (GLsizei)winwidth, (GLsizei)winheight);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(left, left + (height * winwidth) / winheight, bottom, bottom + height, 0, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void PanZoom2D::export_rib(ostream &output) {
  // no projection matrix
  output << "Clipping 1 2000" << endl;     // somewhat arbitrary - hopefully this is plenty of space
  output << "ReverseOrientation" << endl;  // RenderMan has a different handedness from OpenGL's default
  output << "Scale 1 1 -1" << endl;        // so we need to correct for that here
  // scale so that smaller dimension gets scaled to size 2
  scalar scalefactor;
  if (winwidth > winheight)
    scalefactor = 2.0 / height;
  else
    scalefactor = 2.0 / (winwidth * height / winheight);
  output << "Scale " << scalefactor << " " << scalefactor << " 1" << endl;
  // translate so centre of view gets mapped to (0,0,1000)
  output << "Translate " << -(left + 0.5 * winwidth * height / winheight) << " " << -(bottom + 0.5 * height) << " 1000" << endl;
}

void PanZoom2D::display_screen(void) {
  if (action_mode == ZOOM_IN && moved_since_mouse_down) {
    glColor3d(1, 1, 1);
    glBegin(GL_LINE_STRIP);
    glVertex2i(clickx, winheight - clicky);
    glVertex2i(oldmousex, winheight - clicky);
    glVertex2i(oldmousex, winheight - oldmousey);
    glVertex2i(clickx, winheight - oldmousey);
    glVertex2i(clickx, winheight - clicky);
    glEnd();
  }
}

//=================================================================================

StaticText::StaticText(const char *text_) : text(text_) {}

void StaticText::display(int x, int y) {
  dispx = x;
  dispy = y;
  width = glutBitmapLength(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)text) + 1;
  height = 15;
  glColor3d(0.3, 0.3, 0.3);
  glRasterPos2i(x, y - height + 2);
  for (int i = 0; text[i] != 0; ++i) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
  glColor3d(1, 1, 1);
  glRasterPos2i(x + 1, y - height + 3);
  for (int i = 0; text[i] != 0; ++i) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
}

//=================================================================================

Button::Button(const char *text_, int minwidth_) : status(UNINVOLVED), text(text_), minwidth(minwidth_) {}

void Button::display(int x, int y) {
  dispx = x;
  dispy = y;
  int textwidth = glutBitmapLength(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)text);
  if (textwidth < minwidth)
    width = minwidth + 24;
  else
    width = textwidth + 24;
  height = 17;
  if (status == UNINVOLVED) {
    glColor3d(0.7, 0.7, 0.7);
    glBegin(GL_QUADS);
    glVertex2i(x + 1, y - 1);
    glVertex2i(x + width, y - 1);
    glVertex2i(x + width, y - height + 1);
    glVertex2i(x + 1, y - height + 1);
    glEnd();
    glColor3d(0.3, 0.3, 0.3);
    glLineWidth(1);
    glBegin(GL_LINE_STRIP);
    glVertex2i(x, y - 2);
    glVertex2i(x, y - height);
    glVertex2i(x + width - 1, y - height);
    glEnd();
    glColor3d(0.3, 0.3, 0.3);
  } else {
    if (status == SELECTED)
      glColor3d(0.8, 0.8, 0.8);
    else
      glColor3d(1, 1, 1);
    glBegin(GL_QUADS);
    glVertex2i(x, y - 1);
    glVertex2i(x + width, y - 1);
    glVertex2i(x + width, y - height);
    glVertex2i(x, y - height);
    glEnd();
    glColor3d(0, 0, 0);
  }
  glRasterPos2i(x + (width - textwidth) / 2, y - height + 5);
  for (int i = 0; text[i] != 0; ++i) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
}

bool Button::click(int state, int x, int y) {
  if (state == GLUT_DOWN && x > dispx && x <= dispx + width && y < dispy - 2 && y >= dispy - height) {
    status = HIGHLIGHTED;
    glutPostRedisplay();
    return true;
  } else if (state == GLUT_UP && status != UNINVOLVED) {
    status = UNINVOLVED;
    glutPostRedisplay();
    if (x >= dispx && x < dispx + width && y < dispy - 2 && y >= dispy - height) action();
    return true;
  } else
    return false;
}

void Button::drag(int x, int y) {
  // needs to control highlighting (SELECTED vs. HIGHLIGHTED)
  if (status == SELECTED && x >= dispx && x < dispx + width && y < dispy - 2 && y >= dispy - height) {
    status = HIGHLIGHTED;
    glutPostRedisplay();
  } else if (status == HIGHLIGHTED && !(x >= dispx && x < dispx + width && y < dispy - 2 && y >= dispy - height)) {
    status = SELECTED;
    glutPostRedisplay();
  }
}

//=================================================================================

Slider::Slider(const char *text_, int length_, int position_, int justify_)
    : status(UNINVOLVED), text(text_), length(length_), justify(justify_), position(position_) {}

void Slider::display(int x, int y) {
  dispx = x;
  dispy = y;
  width = glutBitmapLength(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)text);
  if (width < justify) width = justify;
  width += 11 + 6 + length + 1;
  height = 15;
  glColor3d(0.3, 0.3, 0.3);
  glRasterPos2i(x, y - height + 2);
  for (int i = 0; text[i] != 0; ++i) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
  glColor3d(1, 1, 1);
  glRasterPos2i(x + 1, y - height + 3);
  for (int i = 0; text[i] != 0; ++i) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
  scrollxmin = x + width - length - 12;
  scrollxmax = x + width;
  scrollymin = y - height + 1;
  scrollymax = y - 2;
  glColor3d(0.3, 0.3, 0.3);
  glLineWidth(1);
  glBegin(GL_LINE_STRIP);
  glVertex2i(scrollxmin, scrollymax - 1);
  glVertex2i(scrollxmin, scrollymin);
  glVertex2i(scrollxmax - 1, scrollymin);
  glVertex2i(scrollxmax - 1, scrollymax - 1);
  glEnd();
  glColor3d(0.7, 0.7, 0.7);
  glBegin(GL_LINE_STRIP);
  glVertex2i(scrollxmin + 1, scrollymax);
  glVertex2i(scrollxmin + 1, scrollymin + 1);
  glVertex2i(scrollxmax, scrollymin + 1);
  glVertex2i(scrollxmax, scrollymax);
  glEnd();
  if (status == UNINVOLVED) {
    glColor3d(0.3, 0.3, 0.3);
    glBegin(GL_LINE_STRIP);
    glVertex2i(scrollxmin + position + 2, scrollymax - 2);
    glVertex2i(scrollxmin + position + 2, scrollymin + 2);
    glVertex2i(scrollxmin + position + 10, scrollymin + 2);
    glEnd();
    glColor3d(0.7, 0.7, 0.7);
    glBegin(GL_QUADS);
    glVertex2i(scrollxmin + position + 3, scrollymin + 3);
    glVertex2i(scrollxmin + position + 11, scrollymin + 3);
    glVertex2i(scrollxmin + position + 11, scrollymax);
    glVertex2i(scrollxmin + position + 3, scrollymax);
    glEnd();
  } else {  // SELECTED
    glColor3d(1, 1, 1);
    glBegin(GL_QUADS);
    glVertex2i(scrollxmin + position + 2, scrollymin + 2);
    glVertex2i(scrollxmin + position + 11, scrollymin + 2);
    glVertex2i(scrollxmin + position + 11, scrollymax);
    glVertex2i(scrollxmin + position + 2, scrollymax);
    glEnd();
  }
}

bool Slider::click(int state, int x, int y) {
  if (state == GLUT_DOWN && x > scrollxmin + position + 2 && x <= scrollxmin + position + 11 && y < scrollymax - 1 && y >= scrollymin + 2) {
    status = SELECTED;
    clickx = x;
    glutPostRedisplay();
    return true;
  } else if (status != UNINVOLVED && state == GLUT_UP) {
    status = UNINVOLVED;
    glutPostRedisplay();
    return true;
  } else
    return false;
}

void Slider::drag(int x, int y) {
  if (status == SELECTED) {
    glutPostRedisplay();
    int newposition = position + (x - clickx);
    clickx = x;
    if (newposition < 0) {
      clickx += (0 - newposition);
      newposition = 0;
    } else if (newposition > length) {
      clickx += (length - newposition);
      newposition = length;
    }
    if (newposition != position) {
      position = newposition;
      action();
      glutPostRedisplay();
    }
  }
}

//=================================================================================

WidgetList::WidgetList(int indent_, bool hidden_) : indent(indent_), hidden(hidden_), downclicked_member(-1) {}

void WidgetList::display(int x, int y) {
  dispx = x;
  dispy = y;
  if (hidden) {
    width = height = 0;
  } else {
    height = 0;
    for (unsigned int i = 0; i < list.size(); ++i) {
      list[i]->display(x + indent, y - height);
      height += list[i]->height;
      width = (width < indent + list[i]->width) ? indent + list[i]->width : width;
    }
  }
}

bool WidgetList::click(int state, int x, int y) {
  // if(hidden || x<dispx || x>=dispx+width || y>=dispy || y<dispy-height)
  // return false; // early exit
  if (state == GLUT_DOWN) {  // search for correct widget
    for (unsigned int i = 0; i < list.size(); ++i) {
      if (list[i]->click(state, x, y)) {
        downclicked_member = i;
        return true;
      }
    }
  } else if (state == GLUT_UP && downclicked_member >= 0) {
    list[downclicked_member]->click(state, x, y);
    downclicked_member = -1;
  }
  return false;
}

void WidgetList::drag(int x, int y) {
  if (downclicked_member >= 0) list[downclicked_member]->drag(x, y);
}

//=================================================================================

static void gluviReshape(int w, int h) {
  winwidth = w;
  winheight = h;
  glutPostRedisplay();  // triggers the camera to adjust itself to the new
                        // dimensions
}

//=================================================================================

static void gluviDisplay() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // draw the scene
  if (camera) camera->gl_transform();
  if (userDisplayFunc) userDisplayFunc();

  // now draw widgets on top
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glLineWidth(1);
  // and probably more needs setting before widgets

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, winwidth, 0, winheight);

  root.display(0, winheight);

  // and allow the camera to draw something on screen (e.g. for zooming extent)
  if (camera) camera->display_screen();

  glutSwapBuffers();
}

static void gluviKeyboard(unsigned char key, int x, int y) {
  if (key == 27 || key == 'q') {
    exit(0);
  } else if (userKeyFunc) {
    userKeyFunc(key, x, y);
  }
}

//=================================================================================

static enum { NOBODY, CAMERA, WIDGETS, USER } mouse_owner = NOBODY;

static void gluviMouse(int button, int state, int x, int y) {
  if (state == GLUT_DOWN) {
    int mods = glutGetModifiers();
    if (camera && mods == GLUT_ACTIVE_SHIFT) {
      camera->click(button, state, x, y);
      mouse_owner = CAMERA;
    } else if (button == GLUT_LEFT_BUTTON && root.click(state, x, winheight - y)) {
      mouse_owner = WIDGETS;
    } else if (userMouseFunc) {
      userMouseFunc(button, state, x, y);
      mouse_owner = USER;
    }
  } else {  // mouse up - send event to whoever got the mouse down
    switch (mouse_owner) {
      case CAMERA:
        camera->click(button, state, x, y);
        break;
      case WIDGETS:
        root.click(state, x, winheight - y);
        break;
      case USER:
        if (userMouseFunc) userMouseFunc(button, state, x, y);
        break;
      default:;  // nothing to do
    }
    mouse_owner = NOBODY;
  }
}

//=================================================================================

static void gluviDrag(int x, int y) {
  switch (mouse_owner) {
    case CAMERA:
      camera->drag(x, y);
      break;
    case WIDGETS:
      root.drag(x, winheight - y);
      break;
    case USER:
      if (userDragFunc) userDragFunc(x, y);
      break;
    default:;  // nothing to do
  }
}

//=================================================================================

void ppm_screenshot(const char *filename_format, ...) {
  va_list ap;
  va_start(ap, filename_format);
#ifdef _MSC_VER
#define FILENAMELENGTH 256
  char filename[FILENAMELENGTH];
  _vsnprintf(filename, FILENAMELENGTH, filename_format, ap);
  ofstream out(filename, ofstream::binary);
#else
  char *filename;
  vasprintf(&filename, filename_format, ap);
  ofstream out(filename, ofstream::binary);
  free(filename);
#endif
  if (!out) return;
  GLubyte *image_buffer = new GLubyte[3 * winwidth * winheight];
  glReadBuffer(GL_FRONT);
  glReadPixels(0, 0, winwidth, winheight, GL_RGB, GL_UNSIGNED_BYTE, image_buffer);
  out << "P6\n" << winwidth << ' ' << winheight << " 255\n";
  for (int i = 1; i <= winheight; ++i) out.write((const char *)image_buffer + 3 * winwidth * (winheight - i), 3 * winwidth);
  delete[] image_buffer;
}

static void write_big_endian_ushort(std::ostream &output, unsigned short v) {
  output.put((v >> 8) % 256);
  output.put(v % 256);
}

static void write_big_endian_uint(std::ostream &output, unsigned int v) {
  output.put((v >> 24) % 256);
  output.put((v >> 16) % 256);
  output.put((v >> 8) % 256);
  output.put(v % 256);
}

void sgi_screenshot(const char *filename_format, ...) {
  va_list ap;
  va_start(ap, filename_format);
#ifdef _MSC_VER
#define FILENAMELENGTH 256
  char filename[FILENAMELENGTH];
  _vsnprintf(filename, FILENAMELENGTH, filename_format, ap);
  ofstream output(filename, ofstream::binary);
#else
  char *filename;
  vasprintf(&filename, filename_format, ap);
  ofstream output(filename, ofstream::binary);
#endif
  if (!output) return;
  // first write the SGI header
  write_big_endian_ushort(output, 474);        // magic number to identify this as an SGI image file
  output.put(0);                               // uncompressed
  output.put(1);                               // use 8-bit colour depth
  write_big_endian_ushort(output, 3);          // number of dimensions
  write_big_endian_ushort(output, winwidth);   // x size
  write_big_endian_ushort(output, winheight);  // y size
  write_big_endian_ushort(output, 3);          // three colour channels (z size)
  write_big_endian_uint(output, 0);            // minimum pixel value
  write_big_endian_uint(output, 255);          // maximum pixel value
  write_big_endian_uint(output, 0);            // dummy spacing
  // image name
  int i;
  for (i = 0; i < 80 && filename[i]; ++i) output.put(filename[i]);
  for (; i < 80; ++i) output.put(0);
  write_big_endian_uint(output, 0);         // colormap is normal
  for (i = 0; i < 404; ++i) output.put(0);  // filler to complete header
  // now write the SGI image data
  GLubyte *image_buffer = new GLubyte[winwidth * winheight];
  glReadBuffer(GL_FRONT);
  glReadPixels(0, 0, winwidth, winheight, GL_RED, GL_UNSIGNED_BYTE, image_buffer);
  output.write((const char *)image_buffer, winwidth * winheight);
  glReadPixels(0, 0, winwidth, winheight, GL_GREEN, GL_UNSIGNED_BYTE, image_buffer);
  output.write((const char *)image_buffer, winwidth * winheight);
  glReadPixels(0, 0, winwidth, winheight, GL_BLUE, GL_UNSIGNED_BYTE, image_buffer);
  output.write((const char *)image_buffer, winwidth * winheight);
  delete[] image_buffer;
#ifndef _MSC_VER
  free(filename);
#endif
}

void set_generic_lights(void) {
  glEnable(GL_LIGHTING);
  {
    float ambient[4] = {.3, .3, .3, 1};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
  }
  {
    float color[4] = {.8, .8, .8, 1};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, color);
    glLightfv(GL_LIGHT0, GL_SPECULAR, color);
    glEnable(GL_LIGHT0);
  }
  {
    float color[4] = {.4, .4, .4, 1};
    glLightfv(GL_LIGHT1, GL_DIFFUSE, color);
    glLightfv(GL_LIGHT1, GL_SPECULAR, color);
    glEnable(GL_LIGHT1);
  }
  {
    float color[4] = {.2, .2, .2, 1};
    glLightfv(GL_LIGHT2, GL_DIFFUSE, color);
    glLightfv(GL_LIGHT2, GL_SPECULAR, color);
    glEnable(GL_LIGHT2);
  }
}

void set_generic_material(scalar r, scalar g, scalar b, GLenum face) {
  float ambient[4], diffuse[4], specular[4];
  ambient[0] = 0.1 * r + 0.03;
  ambient[1] = 0.1 * g + 0.03;
  ambient[2] = 0.1 * b + 0.03;
  ambient[3] = 1;
  diffuse[0] = 0.7 * r;
  diffuse[1] = 0.7 * g;
  diffuse[2] = 0.7 * b;
  diffuse[3] = 1;
  specular[0] = 0.1 * r + 0.1;
  specular[1] = 0.1 * g + 0.1;
  specular[2] = 0.1 * b + 0.1;
  specular[3] = 1;
  glMaterialfv(face, GL_AMBIENT, ambient);
  glMaterialfv(face, GL_DIFFUSE, diffuse);
  glMaterialfv(face, GL_SPECULAR, specular);
  glMaterialf(face, GL_SHININESS, 32);
}

void set_matte_material(scalar r, scalar g, scalar b, GLenum face) {
  float ambient[4], diffuse[4], specular[4];
  ambient[0] = 0.1 * r + 0.03;
  ambient[1] = 0.1 * g + 0.03;
  ambient[2] = 0.1 * b + 0.03;
  ambient[3] = 1;
  diffuse[0] = 0.9 * r;
  diffuse[1] = 0.9 * g;
  diffuse[2] = 0.9 * b;
  diffuse[3] = 1;
  specular[0] = 0;
  specular[1] = 0;
  specular[2] = 0;
  specular[3] = 1;
  glMaterialfv(face, GL_AMBIENT, ambient);
  glMaterialfv(face, GL_DIFFUSE, diffuse);
  glMaterialfv(face, GL_SPECULAR, specular);
}

// draw_2d_arrow assumptions:
// line width, point size, and color are set by the user prior to calling the
// routine
void draw_2d_arrow(const Vector2s base, const Vector2s point, scalar arrow_head_length) {
  Vector2s w = point - base;
  double len = w.norm();

  if (len == 0) {
    glBegin(GL_POINTS);
    glVertex2f(base[0], base[1]);
    glEnd();
    return;
  }

  w = w / (scalar)len;  // normalize to build coordinate system

  // u = w + 90
  // using rotation matrix  0  1
  //	                     -1  0
  Vector2s u = Vector2s(1 * w[1], -1 * w[0]);
  u.normalize();

  // v = w - 90 (in fact v=-u)
  Vector2s v = Vector2s(-1 * w[1], 1 * w[0]);
  v.normalize();

  if (!arrow_head_length) {
    arrow_head_length = 0.1 * len;
  }

  // arrow head points
  Vector2s arrow1, arrow2;
  arrow1 = point + arrow_head_length * (v - w);
  arrow2 = point + arrow_head_length * (u - w);

  glBegin(GL_LINES);
  glVertex2f(base[0], base[1]);
  glVertex2f(point[0], point[1]);
  glVertex2f(point[0], point[1]);
  glVertex2f(arrow1[0], arrow1[1]);
  glVertex2f(point[0], point[1]);
  glVertex2f(arrow2[0], arrow2[1]);
  glEnd();
}

void draw_text(const scalar point[3], const char *text, int fontsize) {
  // please implement me!
}

//=================================================================================

void init(const char *windowtitle, int *argc, char **argv) {
  glutInit(argc, argv);
  glutInitWindowSize(winwidth, winheight);
  glutCreateWindow(windowtitle);
  glutReshapeFunc(gluviReshape);
  glutDisplayFunc(gluviDisplay);
  glutKeyboardFunc(gluviKeyboard);
  glutMouseFunc(gluviMouse);
  glutMotionFunc(gluviDrag);
  glEnable(GL_DEPTH_TEST);
  glClearColor(0, 0, 0, 0);
  glClearDepth(1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
}

//=================================================================================

void (*userDisplayFunc)(void) = 0;
void (*userMouseFunc)(int button, int state, int x, int y) = 0;
void (*userDragFunc)(int x, int y) = 0;
void (*userKeyFunc)(unsigned char key, int x, int y) = 0;
Camera *camera = 0;
WidgetList root(0);
int winwidth = 1536, winheight = 1152;

//=================================================================================

void run(void) { glutMainLoop(); }

};  // namespace Gluvi
