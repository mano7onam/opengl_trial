#define MACOS

#ifdef UBUNTU
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <functional>
#include <array>
#include <cassert>
#include <tuple>
#include <cmath>

#include "kdtree.h"
#include "SetPairs.h"
#include "Triangulation.h"
#include "DataStructure.h"
#include "DSU.h"
#include "math/point.h"
#include "algorithms/NurbsBuilder.h"

using std::vector;
using std::pair;
using std::string;

const bool DRAW_POINTS = false;
const bool DRAW_TRIANGLES = false;
const bool DRAW_NURBS = true;

vector< Vertex > vertsBunny;
std::vector<Triangle> trianglesBunny;
vector<vector<GLfloat>> bunnyNurbs;
GLUnurbsObj *theNurb;

#ifdef UBUNTU
std::string os = "ubuntu";
#else
std::string os = "macOS";
#endif

std::string who = "bunny";

namespace gls {
    float angle = 0;
    float zoom = 0;
    float theta = 0;
    float phi = 0;
    bool gluLookAt_On = true;
    int width;
    int height;
}

void timer( int extra )
{
    // spin
    gls::angle += 0.5;

    glutPostRedisplay();
    glutTimerFunc( 16, timer, 0 );
}

void readVerticesFromPly() {
    std::mt19937::result_type seed = 7147;
    std::mt19937 gen(seed);
    auto real_rand = std::bind(std::uniform_real_distribution<float>(0,1), gen);

#ifdef UBUNTU
    auto in = std::ifstream("/home/mano/CLionProjects/opengl_trial/bunny/bun000.ply");
#else
    vector<string> files = {"../bunny/bun000.ply", /*"../bunny/bun045.ply",*/ "../bunny/bun090.ply", "../bunny/bun180.ply", "../bunny/bun270.ply"/*, "../bunny/bun315.ply"*/};
    vector<float> angs = {0, /*M_PI / 4.0,*/ -M_PI / 2.0, M_PI, -3.0 * M_PI / 2.0};
    vector<int> nums = {40256, /*40097, */30379, 40251, 31701/*, 35336*/};
#endif
    for (int jj = 0; jj < 4/*files.size()*/; ++jj) {
        auto in = std::ifstream(files[jj]);
        std::string dummy;
        for (int i = 0; i < 24; ++i) {
            getline(in, dummy);
        }
        float x, y, z;
        for (int i = 0; i < nums[jj]; ++i) {
            in >> x >> y >> z;
            Vertex cur;
            cur.x = x;
            cur.z = y;
            cur.y = z;

            cur.x = x * cos(angs[jj]) - z * sin(angs[jj]);
            cur.y = x * sin(angs[jj]) + z * cos(angs[jj]);

            cur.a = 1.0f;
            cur.r = 1.0f;
            cur.g = 1.0f;
            cur.b = 1.0f;
            vertsBunny.push_back(cur);
        }
    }

    NurbsBuilder nb(vertsBunny);
    nb.buildNurbs();
    vertsBunny = nb.getNewVertices();
    trianglesBunny = nb.getTriangles();
    bunnyNurbs = nb.getNurbs();
}

template<class T>
void relax_min(T& a, T b) {
    a = std::min(a, b);
}

template<class T>
void relax_max(T &a, T b) {
    a = std::max(a, b);
}

void printMinMaxCoord(std::vector<Vertex> v) {
    float minx = 10000, maxx = -10000, miny = 10000, maxy = -10000, minz = 10000, maxz = -10000;
    float cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < v.size(); ++i) {
        relax_min(minx, v[i].x);
        relax_max(maxx, v[i].x);

        relax_min(miny, v[i].y);
        relax_max(maxy, v[i].y);

        relax_min(minz, v[i].z);
        relax_max(maxz, v[i].z);
    }
    cx = (minx + maxx) / 2.0;
    cy = (miny + maxy) / 2.0;
    cz = (minz + maxz) / 2.0;

    std::cout << minx << ' ' << maxx << '\n';
    std::cout << miny << ' ' << maxy << '\n';
    std::cout << minz << ' ' << maxz << '\n';
    std::cout << cx << ' ' << cy << ' ' << cz << '\n';
    std::cout << '\n';
}

void display(void)
{

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    const double w = glutGet( GLUT_WINDOW_WIDTH );
    const double h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 60.0, w / h, 0.001, 1000.0 );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt( 0.1, 0.3, 0.2, 0, 0, 0.1, 0, 0, 1 );

    glRotatef( gls::angle, 0, 0, 1 );

    if (DRAW_TRIANGLES) {
        for (auto tr : trianglesBunny) {
            glBegin(GL_TRIANGLES);

            glColor3f( tr.a.r, tr.a.g, tr.a.b );
            glVertex3f( tr.a.x, tr.a.y, tr.a.z );

            glColor3f( tr.b.r, tr.b.g, tr.b.b );
            glVertex3f( tr.b.x, tr.b.y, tr.b.z );

            glColor3f( tr.c.r, tr.c.g, tr.c.b );
            glVertex3f( tr.c.x, tr.c.y, tr.c.z );

            glEnd();
        }
    }

    if (DRAW_POINTS) {
        for (auto v : vertsBunny) {
            glBegin(GL_POINTS);
            glColor3f(v.r, v.g, v.b);
            glPointSize(0.000000001);
            glVertex3f(v.x, v.y, v.z);
            glEnd();
        }
    }

    if (DRAW_NURBS) {
        GLfloat knots[12] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        for (auto &nurb : bunnyNurbs) {
            gluBeginSurface(theNurb);
            gluNurbsSurface(theNurb,
                            12, knots, 12, knots,
                            6 * 3, 3, &nurb[0],
                            6, 6, GL_MAP2_VERTEX_3);
            gluEndSurface(theNurb);
        }
    }

    printMinMaxCoord(vertsBunny);

    glutSwapBuffers();
}

// tries to move camera...
void changeSize(int w, int h){
    gls::width = w;
    gls::height = h;

    if(h == 0) h = 1;
    float ratio = 1.0 * w / h;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

//    glViewport(0, 0, w, h);

//    gluPerspective(45, ratio, 1, 1000);

//    float r = 5.0f;
//    float eyeX = r * sin(gls::theta * gls::radianFactor) * cos(gls::phi * gls::radianFactor);
//    float eyeY = r * sin(gls::theta * radianFactor) * sin(gls::phi * radianFactor);
//    float eyeZ = r * cos(radianFactor * gls::theta);
//
//    float centerX = 0, centerY = 0, centerZ = 0;
//    float upX = 0, upY = 1.0f, upZ = 0;
//
//    if(gls::gluLookAt_On) {
//        gluLookAt(eyeX, eyeY, eyeZ,
//                  centerX, centerY, centerZ,
//                  upX, upY, upZ);
//    }

    glScalef(gls::zoom, gls::zoom, gls::zoom);
    glMatrixMode(GL_MODELVIEW);
}

void inputKey(unsigned char c, int x, int y){
    switch (c) {
        case 'i' : gls::zoom = gls::zoom+ 0.5; break;
        case 'o' : gls::zoom = gls::zoom-0.5; break;
        case 't' : gls::theta++; if(gls::theta > 360) gls::theta = 1; break;
        case 'p' : gls::phi++; if(gls::phi > 360) gls::phi = 1; break;
        case 'T' : gls::theta--; if(gls::theta < 0) gls::theta = 359; break;
        case 'P' : gls::phi--; if(gls::phi < 0) gls::phi = 359; break;
        case 'g' : gls::gluLookAt_On = !gls::gluLookAt_On;; break;
    }
    changeSize(gls::width, gls::height);
}

void nurbsError(GLenum errorCode)
{
    const GLubyte *estring;

    estring = gluErrorString(errorCode);
    fprintf (stderr, "Nurbs Error: %s\n", estring);
    exit (0);
}

int main( int argc, char **argv )
{
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );
    glutInitWindowSize( 800, 600 );
    glutCreateWindow( "Attractor" );

    glutDisplayFunc( display );
    glutKeyboardFunc(inputKey);
    glutReshapeFunc(changeSize);
    glutTimerFunc( 0, timer, 0 );

    GLfloat mat_diffuse[] = { 0.7, 0.7, 0.7, 1.0 };
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 100.0 };

    glClearColor (0.0, 0.0, 0.0, 0.0);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    theNurb = gluNewNurbsRenderer();
    gluNurbsProperty(theNurb, GLU_SAMPLING_TOLERANCE, 25.0);
    gluNurbsProperty(theNurb, GLU_DISPLAY_MODE, GLU_FILL);
    gluNurbsCallback(theNurb, GLU_ERROR,
                     (GLvoid (*)()) nurbsError);

    readVerticesFromPly();

    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    glEnable( GL_POINT_SMOOTH );
    glPointSize(1.0f);

    glutMainLoop();
    return 0;
}
