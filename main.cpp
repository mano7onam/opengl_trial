#define UBUNTU

#ifdef UBUNTU
#include <GL/glut.h>
#include "Triangulation.h"
#include "DataStructure.h"
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

#include "kdtree.h"
#include "SetPairs.h"

using std::vector;
using std::pair;

const int K = 50;

struct Vertex
{
    float x, y, z, w;
    float r, g, b, a;
};
vector< Vertex > verts;
vector< Vertex > vertsBunny;
vector<std::vector<int>> bunnyG;
vector<std::tuple<int, int, int>*> bunnyMesh;

class MyPoint : public std::array<float, 3>
{
public:

    // dimension of space (or "k" of k-d tree)
    // KDTree class accesses this member
    static const int DIM = 3;

    // the constructors
    MyPoint() {}
    MyPoint(float x, float y, float z)
    {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
    }
};

vector<MyPoint> myPointsBunny;
kdt::KDTree<MyPoint> kdTree;
SetPairs borderEdges;

void buildKdTreeFromVerts() {
    for (auto v : vertsBunny) {
        myPointsBunny.emplace_back(v.x, v.y, v.z);
    }
    kdTree.build(myPointsBunny);
}

#define sqr(x) ((x) * (x))

float getDist(Vertex a, Vertex b) {
    return sqrt(sqr(a.x - b.x) + sqr(a.y - b.y) + sqr(a.z - b.z));
}

float getDist(MyPoint a, MyPoint b) {
    return sqrt(sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2]));
}

MyPoint getMid(MyPoint a, MyPoint b) {
    return MyPoint((a[0] + b[0]) * 0.5f, (a[1] + b[1]) * 0.5f, (a[2] + b[2]) * 0.5f);
}

struct Triangle {
    Vertex a, b, c;

    Triangle(Vertex a, Vertex b, Vertex c): a(a), b(b), c(c) {}
};
std::vector<Triangle> trianglesBunny;

#ifdef UBUNTU
std::string os = "ubuntu";
#else
std::string os = "macOS";
#endif

std::string who = "bunny";

void fillVerts()
{
    // calculate vertices
    // http://paulbourke.net/fractals/lorenz/
    double h = 0.01;
    double a = 10.0;
    double b = 28.0;
    double c = 8.0 / 3.0;

    Vertex cur;
    cur.a = 0.09f;

    double x0 = 0.1;
    double y0 = 0;
    double z0 = 0;
    for( unsigned int i = 0; i < 100000; i++ )
    {
        if(i == 20000)
        {
            cur.r = 1.0f;
            cur.g = 0.0f;
            cur.b = 0.0f;
        }
        if(i == 40000)
        {
            cur.r = 1.0f;
            cur.g = 0.0f;
            cur.b = 1.0f;
        }
        if(i == 60000)
        {
            cur.r = 0.0f;
            cur.g = 0.0f;
            cur.b = 1.0f;
        }
        if(i == 80000)
        {
            cur.r = 0.0f;
            cur.g = 1.0f;
            cur.b = 1.0f;
        }

        const double x1 = x0 + h * a * (y0 - x0);
        const double y1 = y0 + h * (x0 * (b - z0) - y0);
        const double z1 = z0 + h * (x0 * y0 - c * z0);
        x0 = x1;
        y0 = y1;
        z0 = z1;

        if( i > 100 )
        {
            cur.x = x0;
            cur.y = y0;
            cur.z = z0;
            verts.push_back( cur );
        }
    }
}

float angle = 0;
void timer( int extra )
{
    // spin
    angle += 0.5;

    glutPostRedisplay();
    glutTimerFunc( 16, timer, 0 );
}

void buildTrianglesBunnyOld() {
    buildKdTreeFromVerts();

    bunnyG.assign(vertsBunny.size(), std::vector<int>());
    for (int i = 0; i < vertsBunny.size(); ++i) {
        std::cerr << i << std::endl;
//        vector<pair<float, int>> ps;
//        for (int j = 0; j < vertsBunny.size(); ++j) {
//            if (i == j) {
//                continue;
//            }
//            ps.emplace_back(getDist(vertsBunny[i], vertsBunny[j]), j);
//        }
//        std::sort(ps.begin(), ps.end());
//        for (int j = 0; j < std::min(K, static_cast<int>(ps.size())); ++j) {
//            bunnyG[i].push_back(ps[j].second);
//        }
        bunnyG[i] = kdTree.knnSearch(myPointsBunny[i], K);
    }

    for (int i = 0; i < bunnyG.size(); ++i) {
        for (int j = 0; j < bunnyG[i].size(); ++j) {
            int n1 = bunnyG[i][j];
            for (int k = 0; k < bunnyG[n1].size(); ++k) {
                int n2 = bunnyG[n1][k];
                if (n2 == i || n2 == n1) {
                    continue;
                }
                if (std::find(bunnyG[i].begin(), bunnyG[i].end(), n2) != bunnyG[i].end()) {
                    trianglesBunny.emplace_back(vertsBunny[i], vertsBunny[n1], vertsBunny[n2]);
                }
            }
        }
    }
}

int nNeibhors[] = {3, 4, 5, 10, 20, 50};
void connectEdgeWithNewVert(pair<int, int> p) {
    auto a = myPointsBunny[p.first];
    auto b = myPointsBunny[p.second];
    auto mab = getMid(a, b);
    int found = -1;
    for (int nb : nNeibhors) {
        vector<int> inds = kdTree.knnSearch(mab, nb);
        for (int id : inds) {
            if (id == p.first || id == p.second) continue;
            if (borderEdges.havePair({p.first, id}) && borderEdges.havePair({p.second, id})) continue;
            if (borderEdges.havePair(found) && !borderEdges.havePair({p.first, id}) && !borderEdges.havePair({p.second, id})) continue;
            found = id;
            break;
        }
        if (found != -1) {
            break;
        }
    }
    if (found != -1) {
        if (borderEdges.havePair({p.first, found})) {
            borderEdges.erasePair({p.first, found});
            borderEdges.addPair({p.second, found});
        } else if (borderEdges.havePair({p.second, found})) {
            borderEdges.erasePair({p.second, found});
            borderEdges.addPair({p.first, found});
        } else {
            borderEdges.addPair({p.second, found});
            borderEdges.addPair({p.first, found});
        }
        trianglesBunny.emplace_back(vertsBunny[p.first], vertsBunny[p.second], vertsBunny[found]);
    }
    borderEdges.erasePair(p);
}

void buildTrianglesBunnyOld1() {
    assert(vertsBunny.size() >= 3);

    buildKdTreeFromVerts();

    std::mt19937::result_type seed = 7147;
    std::mt19937 gen(seed);
    auto indexRand = std::bind(std::uniform_int_distribution<int>(0, static_cast<int>(vertsBunny.size()) - 1), gen);

    int starti = indexRand();
    auto vnbi = kdTree.knnSearch(myPointsBunny[starti], 2);
    int nbi = vnbi[0];
    if (vnbi[0] == starti) nbi = vnbi[1];

    borderEdges.addPair({starti, nbi});
    int counter = 0;
    while (!borderEdges.isEmpty()/* && counter < 10000*/) {
        counter++;
        auto p = borderEdges.getFirstPair();
        connectEdgeWithNewVert(p);
    }
}

void buildTrianglesBunny() {
    vector<dt::Vector3D*> dots;
    for (auto v : vertsBunny) {
        dots.push_back(new dt::Vector3D(v.x, v.y, v.z));
    }
    dt::DelaunayTriangulation delanuaryTriangulation;
    bunnyMesh = delanuaryTriangulation.GetTriangulationResult(dots);
}

void readVerticesFromPly() {
    std::mt19937::result_type seed = 7147;
    std::mt19937 gen(seed);
    auto real_rand = std::bind(std::uniform_real_distribution<float>(0,1), gen);

#ifdef UBUNTU
    auto in = std::ifstream("/home/mano/CLionProjects/opengl_trial/bunny/bun000.ply");
#else
    auto in = std::ifstream("../bunny/bun000.ply");
#endif
    std::string dummy;
    for (int i = 0; i < 24; ++i) {
        getline(in, dummy);
    }
    float x, y, z;
    for (int i = 0; i < 40256; ++i) {
        in >> x >> y >> z;
        Vertex cur;
        cur.x = x;
        cur.y = z;
        cur.z = y;
        cur.a = 0.5f;
        cur.r = real_rand();
        cur.g = real_rand();
        cur.b = real_rand();
//        cur.r = 0.5;
//        cur.g = 0.9;
//        cur.b = 0.2;
        vertsBunny.push_back(cur);
    }

    buildTrianglesBunny();
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
//    gluPerspective( 60.0, w / h, 0.001, 1.0 );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
//    gluLookAt( 70, 70, 70, 0, 0, 0, 0, 0, 1 );
    gluLookAt( 0.1, 0.3, 0.2, 0, 0, 0.1, 0, 0, 1 );

    glRotatef( angle, 0, 0, 1 );

    // draw curve
//    glEnableClientState( GL_VERTEX_ARRAY );
//    glEnableClientState( GL_COLOR_ARRAY );
//    if (who == std::string("bunny")) {
//        glVertexPointer( 3, GL_FLOAT, sizeof( Vertex ), &verts_bunny[0].x );
//        glColorPointer( 4, GL_FLOAT, sizeof( Vertex ), &verts_bunny[0].r );
//        glDrawArrays( GL_LINE_STRIP, 0, verts_bunny.size() );
//    } else {
//        glVertexPointer( 3, GL_FLOAT, sizeof( Vertex ), &verts[0].x );
//        glColorPointer( 4, GL_FLOAT, sizeof( Vertex ), &verts[0].r );
//        glDrawArrays( GL_LINE_STRIP, 0, verts.size() );
//    }
//    glDisableClientState( GL_VERTEX_ARRAY );
//    glDisableClientState( GL_COLOR_ARRAY );

//    for (auto tr : trianglesBunny) {
//        glBegin(GL_TRIANGLES);
//
//        glColor3f( tr.a.r, tr.a.g, tr.a.b );
//        glVertex3f( tr.a.x, tr.a.y, tr.a.z );
//
//        glColor3f( tr.b.r, tr.b.g, tr.b.b );
//        glVertex3f( tr.b.x, tr.b.y, tr.b.z );
//
//        glColor3f( tr.c.r, tr.c.g, tr.c.b );
//        glVertex3f( tr.c.x, tr.c.y, tr.c.z );
//
//        glEnd();
//    }

    for (auto tr : bunnyMesh) {
        glBegin(GL_TRIANGLES);

        Vertex a = vertsBunny[std::get<0>(*tr)];
        Vertex b = vertsBunny[std::get<1>(*tr)];
        Vertex c = vertsBunny[std::get<2>(*tr)];

        glColor3f( a.r, a.g, a.b );
        glVertex3f( a.x, a.y, a.z );

        glColor3f( b.r, b.g, b.b );
        glVertex3f( b.x, b.y, b.z );

        glColor3f( c.r, c.g, c.b );
        glVertex3f( c.x, c.y, c.z );

        glEnd();
    }

    printMinMaxCoord(verts);
    printMinMaxCoord(vertsBunny);

    glutSwapBuffers();
}

int main( int argc, char **argv )
{
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );
    glutInitWindowSize( 800, 600 );
    glutCreateWindow( "Attractor" );

    glutDisplayFunc( display );
    glutTimerFunc( 0, timer, 0 );

    fillVerts();
    readVerticesFromPly();

    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    glEnable( GL_POINT_SMOOTH );
    glPointSize(1.0f);

    glutMainLoop();
    return 0;
}
