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

const int K = 50;

vector< Vertex > verts;
vector< Vertex > vertsBunny;
vector<int> vertsBunnyClasterId;
vector<vector<int>> clastersVertsBunny; // to draw different colors
vector<std::vector<int>> bunnyG;
DSU bunnyDSU;
vector<std::tuple<int, int, int>*> bunnyMesh;
vector<vector<GLfloat>> bunnyNurbs;
GLUnurbsObj *theNurb;

vector<MyPoint> myPointsBunny;
kdt::KDTree<MyPoint> kdTree;
SetPairs bunnyBorderEdges;

void buildKdTreeFromVerts() {
    for (auto v : vertsBunny) {
        myPointsBunny.emplace_back(v.x, v.y, v.z);
    }
    kdTree.build(myPointsBunny);
}


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

//int nNeighbors[] = {3, 4, 5, 10, 20, 50, 100};
//void connectEdgeWithNewVertOld(pair<int, int> p) {
//    auto a = myPointsBunny[p.first];
//    auto b = myPointsBunny[p.second];
//    auto mab = getMid(a, b);
//    int found = -1;
//    for (int nb : nNeighbors) {
//        vector<int> inds = kdTree.knnSearch(mab, nb);
//        for (int id : inds) {
//            if (id == p.first || id == p.second) continue;
//            if (bunnyBorderEdges.havePair({p.first, id}) && bunnyBorderEdges.havePair({p.second, id})) continue;
//            if (bunnyBorderEdges.havePair(found) && !bunnyBorderEdges.havePair({p.first, id}) && !bunnyBorderEdges.havePair({p.second, id})) continue;
//            found = id;
//            break;
//        }
//        if (found != -1) {
//            break;
//        }
//    }
//    if (found != -1) {
//        if (bunnyBorderEdges.havePair({p.first, found})) {
//            bunnyBorderEdges.erasePair({p.first, found});
//            bunnyBorderEdges.addPair({p.second, found});
//        } else if (bunnyBorderEdges.havePair({p.second, found})) {
//            bunnyBorderEdges.erasePair({p.second, found});
//            bunnyBorderEdges.addPair({p.first, found});
//        } else {
//            bunnyBorderEdges.addPair({p.second, found});
//            bunnyBorderEdges.addPair({p.first, found});
//        }
//        trianglesBunny.emplace_back(vertsBunny[p.first], vertsBunny[p.second], vertsBunny[found]);
//    }
//    bunnyBorderEdges.erasePair(p);
//}
//
//void buildTrianglesBunnyOld1() {
//    assert(vertsBunny.size() >= 3);
//
//    buildKdTreeFromVerts();
//
//    std::mt19937::result_type seed = 7147;
//    std::mt19937 gen(seed);
//    auto indexRand = std::bind(std::uniform_int_distribution<int>(0, static_cast<int>(vertsBunny.size()) - 1), gen);
//
//    int starti = indexRand();
//    auto vnbi = kdTree.knnSearch(myPointsBunny[starti], 2);
//    int nbi = vnbi[0];
//    if (vnbi[0] == starti) nbi = vnbi[1];
//
//    bunnyBorderEdges.addPair({starti, nbi});
//    int counter = 0;
//    while (!bunnyBorderEdges.isEmpty()/* && counter < 10000*/) {
//        counter++;
//        auto p = bunnyBorderEdges.getFirstPair();
//        connectEdgeWithNewVert(p);
//    }
//}

void buildTrianglesBunnyDelaunay() {
    vector<dt::Vector3D*> dots;
    for (auto v : vertsBunny) {
        dots.push_back(new dt::Vector3D(v.x, v.y, v.z));
    }
    dt::DelaunayTriangulation delaunayTriangulation;
    bunnyMesh = delaunayTriangulation.GetTriangulationResult(dots);
    int aaaaa = 145;
}

vector<int> getNotUsedTriangle(vector<int> &nbs, vector<bool> &used) {
    vector<int> result;
    int cur = 0;
    for (int nb : nbs) {
        if (used[nb]) {
            continue;
        }
        result.push_back(nb);
        if (result.size() == 3U) {
            break;
        }
    }
    return result;
}

bool haveNeighbor(int a, int b) {
    for (int v : bunnyG[a]) {
        if (v == b) {
            return true;
        }
    }
    return false;
}

void connectVertices(int a, int b) {
    bunnyG[a].push_back(b);
    bunnyG[b].push_back(a);
}

int nNeighbors[] = {3, 4, 5, 6, 7, 8, 10};
void connectEdgeWithNewVert(pair<int, int> p) {
    auto a = myPointsBunny[p.first];
    auto b = myPointsBunny[p.second];
    auto mab = getMid(a, b);
    int found = -1;

    for (int nb : nNeighbors) {
        vector<int> inds = kdTree.knnSearch(mab, nb);
        for (int id : inds) {
            if (id == p.first || id == p.second) continue;
            if (haveNeighbor(p.first, id) && haveNeighbor(p.second, id)) continue;
//            if (bunnyG[id].size() > 6) continue;

            int pid = bunnyDSU.findSet(id);
            int pa = bunnyDSU.findSet(p.first);
            int pb = bunnyDSU.findSet(p.second);
            assert(pa == pb);
            if (pid == pa && !haveNeighbor(id, p.first) && !haveNeighbor(id, p.second)) continue;

            found = id;
            break;
        }
        if (found != -1) {
            break;
        }
    }

    if (found != -1) {
        if (haveNeighbor(found, p.first)) {
            bunnyBorderEdges.erasePair({p.first, found});
            bunnyBorderEdges.addPair({p.second, found});
            connectVertices(p.second, found);
        } else if (haveNeighbor(found, p.second)) {
            bunnyBorderEdges.erasePair({p.second, found});
            bunnyBorderEdges.addPair({p.first, found});
            connectVertices(p.first, found);
        } else {
            bunnyBorderEdges.addPair({p.second, found});
            bunnyBorderEdges.addPair({p.first, found});
            connectVertices(p.first, found);
            connectVertices(p.second, found);
        }
        trianglesBunny.emplace_back(vertsBunny[p.first], vertsBunny[p.second], vertsBunny[found]);
        bunnyDSU.unionSets(p.first, found);
    }

    bunnyBorderEdges.erasePair(p);
}

void saveTrianglesToFile() {
    std::ofstream out("triangles000.txt");
    out << trianglesBunny.size() << '\n';
    for (int i = 0; i < trianglesBunny.size(); ++i) {
        auto tr = trianglesBunny[i];
        out << tr.a.x << ' ' << tr.a.y << ' ' << tr.a.z << '\n';
        out << tr.b.x << ' ' << tr.b.y << ' ' << tr.b.z << '\n';
        out << tr.c.x << ' ' << tr.c.y << ' ' << tr.c.z << '\n';
    }
}

void readTiranglesFromFile() {
    std::ifstream in("triangles000.txt");
    int n;

}

void buildTrianglesBunny() {
    buildKdTreeFromVerts();

    bunnyDSU.init(vertsBunny.size());
    bunnyG.assign(vertsBunny.size(), vector<int>());

    vector<bool> used(vertsBunny.size(), false);
    for (int i = 0; i < vertsBunny.size(); ++i) {
        if (used[i]) {
            continue;
        }
        for (int nnb : nNeighbors) {
            auto nbs = kdTree.knnSearch(myPointsBunny[i], nnb);
            auto tr = getNotUsedTriangle(nbs, used);
            if (tr.size() == 3U) {
                trianglesBunny.emplace_back(vertsBunny[tr[0]], vertsBunny[tr[1]], vertsBunny[tr[2]]);
                for (int j = 0; j < 3; ++j) {
                    used[tr[j]] = true;
                    for (int k = 0; k < 3; ++k) {
                        if (k == j) {
                            continue;
                        }
                        if (j < k) {
                            bunnyBorderEdges.addPair({tr[j], tr[k]});
                        }
                        bunnyG[tr[j]].push_back(tr[k]);
                        bunnyDSU.unionSets(tr[j], tr[k]);
                    }
                }
                break;
            }
        }
    }

    int counter = 0;
    while (!bunnyBorderEdges.isEmpty()) {
        counter++;
        auto p = bunnyBorderEdges.getFirstPair();
        connectEdgeWithNewVert(p);
    }
}

void colorClasters() {
    std::mt19937::result_type seed = 7147;
    std::mt19937 gen(seed);
    auto colorRand = std::bind(std::uniform_real_distribution<float>(0.0f, 0.97f), gen);

    for (const auto &cl: clastersVertsBunny) {
        float r = colorRand();
        float g = colorRand();
        float b = colorRand();
//        float r = 1.0f;
//        float g = 0.0f;
//        float b = 1.0f;
        for (int ind: cl) {
            auto &v = vertsBunny[ind];
            v.r = r;
            v.g = g;
            v.b = b;
        }
    }
}

void clasterizePoints() {
    buildKdTreeFromVerts();

    vertsBunnyClasterId.assign(vertsBunny.size(), -1);

    const int wantNeib = 200;
    // TODO: find nearest not used points to current (need to find 36 this points) (maybe shuffle points at first)
    for (int i = 0; i < vertsBunny.size(); ++i) {
        if (vertsBunnyClasterId[i] != -1) continue;
        auto nbs = kdTree.knnSearch(myPointsBunny[i], wantNeib);
        vector<int> curClaster;
        for (int j = 0; j < nbs.size(); ++j) {
            int id = nbs[j];
            if (vertsBunnyClasterId[i] == -1) {
                curClaster.push_back(id);
            }
        }
        int rest = wantNeib - (int)curClaster.size();
        if (rest > wantNeib / 10) {
            continue;
        }
        for (int j = 0; j < curClaster.size(); ++j) {
            vertsBunnyClasterId[curClaster[j]] = (int)clastersVertsBunny.size();
        }
        clastersVertsBunny.push_back(curClaster);
    }
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

//            cur.z = y * cos(angs[jj]) - x * sin(angs[jj]);
//            cur.x = y * sin(angs[jj]) + x * cos(angs[jj]);

//            cur.y = cur.y * cos(angs[jj]) - cur.z * sin(angs[jj]);
//            cur.z = cur.y * sin(angs[jj]) + cur.z * cos(angs[jj]);

            cur.a = 1.0f;
            cur.r = 1.0f;
            cur.g = 1.0f;
            cur.b = 1.0f;
            vertsBunny.push_back(cur);
        }
    }

//    buildTrianglesBunny();
//    saveTrianglesToFile();

//    clasterizePoints();
//    colorClasters();
    NurbsBuilder nb(vertsBunny);
    nb.buildNurbs();
//    trianglesBunny = nb.getTriangles();
    vertsBunny = nb.getNewVertices();
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
//    gluPerspective( 60.0, w / h, 0.001, 1.0 );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
//    gluLookAt( 70, 70, 70, 0, 0, 0, 0, 0, 1 );
    gluLookAt( 0.1, 0.3, 0.2, 0, 0, 0.1, 0, 0, 1 );

    glRotatef( gls::angle, 0, 0, 1 );

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

//    for (auto v : vertsBunny) {
//        glBegin(GL_POINTS);
//        glColor3f(v.r,v.g,v.b);
//        glPointSize(0.000000001);
//        glVertex3f(v.x, v.y, v.z);
//        glEnd();
//    }

    float minS = 10000000;
    for (auto tr: bunnyMesh) {
        Vertex a = vertsBunny[std::get<0>(*tr)];
        Vertex b = vertsBunny[std::get<1>(*tr)];
        Vertex c = vertsBunny[std::get<2>(*tr)];
        float s = getS(a, b, c);
        if (s > 0.0) {
            minS = std::min(minS, s);
        }
    }

    GLfloat knots[12] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    for (auto &nurb : bunnyNurbs) {
        gluBeginSurface(theNurb);
        gluNurbsSurface(theNurb,
                        12, knots, 12, knots,
                        6 * 3, 3, &nurb[0],
                        6, 6, GL_MAP2_VERTEX_3);
        gluEndSurface(theNurb);
    }

    for (auto tr : bunnyMesh) {
        Vertex a = vertsBunny[std::get<0>(*tr)];
        Vertex b = vertsBunny[std::get<1>(*tr)];
        Vertex c = vertsBunny[std::get<2>(*tr)];

        float s = getS(a, b, c);

//        if (s > minS * 20) {
//            continue;
//        }

        glBegin(GL_TRIANGLES);

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

    fillVerts();
    readVerticesFromPly();

    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    glEnable( GL_POINT_SMOOTH );
    glPointSize(1.0f);

    glutMainLoop();
    return 0;
}
