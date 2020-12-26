/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2020 by Denis Haumont

This file is part of Visilib.

Visilib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Visilib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Visilib. If not, see <http://www.gnu.org/licenses/>
*/

#include <math.h>
#include <string.h>

static double _left = 0.0;
static double _right = 0.0;
static double _bottom = 0.0;
static double _top = 0.0;
static double _zNear = -10.0;
static double _zFar = 10.0;
static double _scaling = 1.0;
static int  _mouseX = 0;
static int  _mouseY = 0;
static bool _mouseLeft = false;
static bool _mouseMiddle = false;
static bool _mouseRight = false;
static double _dragPosX = 0.0;
static double _dragPosY = 0.0;
static double _dragPosZ = 0.0;

static double vlen(double x, double y, double z);
#ifdef USE_GLUT
GLfloat zprReferencePoint[4] = { 0, 0, 0, 0 };
#endif
static double _matrix[16];
static double _matrixInverse[16];

static void setViewPortScaling(double s);
static void   pos(double *px, double *py, double *pz, int x, int y, int *viewport);

#ifdef USE_GLUT
static void
invertMatrix(const GLdouble *m, GLdouble *out)
{
    /* NB. OpenGL Matrices are COLUMN major. */
#define MAT(m,r,c) (m)[(c)*4+(r)]

    /* Here's some shorthand converting standard (row,column) to index. */
#define m11 MAT(m,0,0)
#define m12 MAT(m,0,1)
#define m13 MAT(m,0,2)
#define m14 MAT(m,0,3)
#define m21 MAT(m,1,0)
#define m22 MAT(m,1,1)
#define m23 MAT(m,1,2)
#define m24 MAT(m,1,3)
#define m31 MAT(m,2,0)
#define m32 MAT(m,2,1)
#define m33 MAT(m,2,2)
#define m34 MAT(m,2,3)
#define m41 MAT(m,3,0)
#define m42 MAT(m,3,1)
#define m43 MAT(m,3,2)
#define m44 MAT(m,3,3)

    GLdouble det;
    GLdouble d12, d13, d23, d24, d34, d41;
    GLdouble tmp[16]; /* Allow out == in. */

    /* Inverse = adjoint / det. (See linear algebra texts.)*/

    /* pre-compute 2x2 dets for last two rows when computing */
    /* cofactors of first two rows. */
    d12 = (m31*m42 - m41*m32);
    d13 = (m31*m43 - m41*m33);
    d23 = (m32*m43 - m42*m33);
    d24 = (m32*m44 - m42*m34);
    d34 = (m33*m44 - m43*m34);
    d41 = (m34*m41 - m44*m31);

    tmp[0] = (m22 * d34 - m23 * d24 + m24 * d23);
    tmp[1] = -(m21 * d34 + m23 * d41 + m24 * d13);
    tmp[2] = (m21 * d24 + m22 * d41 + m24 * d12);
    tmp[3] = -(m21 * d23 - m22 * d13 + m23 * d12);

    /* Compute determinant as early as possible using these cofactors. */
    det = m11 * tmp[0] + m12 * tmp[1] + m13 * tmp[2] + m14 * tmp[3];

    /* Run singularity test. */
    if (det == 0.0) {
        /* printf("invert_matrix: Warning: Singular matrix.\n"); */
        /*    memcpy(out,_identity,16*sizeof(double)); */
    }
    else {
        GLdouble invDet = 1.0 / det;
        /* Compute rest of inverse. */
        tmp[0] *= invDet;
        tmp[1] *= invDet;
        tmp[2] *= invDet;
        tmp[3] *= invDet;

        tmp[4] = -(m12 * d34 - m13 * d24 + m14 * d23) * invDet;
        tmp[5] = (m11 * d34 + m13 * d41 + m14 * d13) * invDet;
        tmp[6] = -(m11 * d24 + m12 * d41 + m14 * d12) * invDet;
        tmp[7] = (m11 * d23 - m12 * d13 + m13 * d12) * invDet;

        /* Pre-compute 2x2 dets for first two rows when computing */
        /* cofactors of last two rows. */
        d12 = m11*m22 - m21*m12;
        d13 = m11*m23 - m21*m13;
        d23 = m12*m23 - m22*m13;
        d24 = m12*m24 - m22*m14;
        d34 = m13*m24 - m23*m14;
        d41 = m14*m21 - m24*m11;

        tmp[8] = (m42 * d34 - m43 * d24 + m44 * d23) * invDet;
        tmp[9] = -(m41 * d34 + m43 * d41 + m44 * d13) * invDet;
        tmp[10] = (m41 * d24 + m42 * d41 + m44 * d12) * invDet;
        tmp[11] = -(m41 * d23 - m42 * d13 + m43 * d12) * invDet;
        tmp[12] = -(m32 * d34 - m33 * d24 + m34 * d23) * invDet;
        tmp[13] = (m31 * d34 + m33 * d41 + m34 * d13) * invDet;
        tmp[14] = -(m31 * d24 + m32 * d41 + m34 * d12) * invDet;
        tmp[15] = (m31 * d23 - m32 * d13 + m33 * d12) * invDet;

        memcpy(out, tmp, 16 * sizeof(GLdouble));
    }

#undef m11
#undef m12
#undef m13
#undef m14
#undef m21
#undef m22
#undef m23
#undef m24
#undef m31
#undef m32
#undef m33
#undef m34
#undef m41
#undef m42
#undef m43
#undef m44
#undef MAT
}

static void
getMatrix()
{
    glGetDoublev(GL_MODELVIEW_MATRIX, _matrix);
    invertMatrix(_matrix, _matrixInverse);
}

static double
vlen(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

static void
pos(double *px, double *py, double *pz, int x, int y, int *viewport)
{
    /*
      Use the ortho projection and viewport information
      to map from mouse co-ordinates back into world
      co-ordinates
      */

    *px = (double)(x - viewport[0]) / (double)(viewport[2]);
    *py = (double)(y - viewport[1]) / (double)(viewport[3]);

    *px = _left + (*px)*(_right - _left);
    *py = _top + (*py)*(_bottom - _top);
    *pz = _zNear;
}

static void
zprMouse(int button, int state, int x, int y)
{
    GLint viewport[4];

    _mouseX = x;
    _mouseY = y;

    if (state == GLUT_UP)
        switch (button)
    {
        case GLUT_LEFT_BUTTON:   _mouseLeft = false; break;
        case GLUT_MIDDLE_BUTTON: _mouseMiddle = false; break;
        case GLUT_RIGHT_BUTTON:  _mouseRight = false; break;
    }
    else
        switch (button)
    {
        case GLUT_LEFT_BUTTON:   _mouseLeft = true; break;
        case GLUT_MIDDLE_BUTTON: _mouseMiddle = true; break;
        case GLUT_RIGHT_BUTTON:  _mouseRight = true; break;
    }

    glGetIntegerv(GL_VIEWPORT, viewport);
    pos(&_dragPosX, &_dragPosY, &_dragPosZ, x, y, viewport);
    glutPostRedisplay();
}

static void
zprMotion(int x, int y)
{
    bool changed = false;

    int dx = x - _mouseX;
    int dy = y - _mouseY;

    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    if (dx == 0 && dy == 0)
        return;

    if (_mouseMiddle || (_mouseLeft && _mouseRight))
    {
        double s = exp((double)dy*0.01);

        glTranslatef(zprReferencePoint[0], zprReferencePoint[1], zprReferencePoint[2]);
        glScaled(s, s, s);
        glTranslated(-zprReferencePoint[0], -zprReferencePoint[1], -zprReferencePoint[2]);

        changed = true;
    }
    else
        if (_mouseLeft)
        {
            double ax, ay, az;
            double bx, by, bz;
            double angle;

            ax = dy;
            ay = dx;
            az = 0.0;
            angle = vlen(ax, ay, az) / (double)(viewport[2] + 1)*180.0;

            /* Use inverse matrix to determine local axis of rotation */

            bx = _matrixInverse[0] * ax + _matrixInverse[4] * ay + _matrixInverse[8] * az;
            by = _matrixInverse[1] * ax + _matrixInverse[5] * ay + _matrixInverse[9] * az;
            bz = _matrixInverse[2] * ax + _matrixInverse[6] * ay + _matrixInverse[10] * az;

            glTranslatef(zprReferencePoint[0], zprReferencePoint[1], zprReferencePoint[2]);
            glRotated(angle, bx, by, bz);
            glTranslatef(-zprReferencePoint[0], -zprReferencePoint[1], -zprReferencePoint[2]);

            changed = true;
        }
        else
            if (_mouseRight)
            {
                double px, py, pz;

                pos(&px, &py, &pz, x, y, viewport);

                glLoadIdentity();
                glTranslated(px - _dragPosX, py - _dragPosY, pz - _dragPosZ);
                glMultMatrixd(_matrix);

                _dragPosX = px;
                _dragPosY = py;
                _dragPosZ = pz;

                changed = true;
            }

    _mouseX = x;
    _mouseY = y;

    if (changed)
    {
        getMatrix();
        glutPostRedisplay();
    }
}

static void
zprReshape(int w, int h)
{
    glViewport(0, 0, w, h);

    _top = _scaling;
    _bottom = -_scaling;

    _left = -_scaling;
    _right = _scaling;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(_left, _right, _bottom, _top, -100*_scaling, 100*_scaling);

    glMatrixMode(GL_MODELVIEW);
}
#endif

static void setViewPortScaling(double s)
{
    _scaling = s;
 #ifdef USE_GLUT
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    zprReshape(viewport[2], viewport[3]);
#endif
}

