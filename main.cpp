#include <GLUT/glut.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <functional>

using std::vector;
using std::pair;

const int K = 5;

struct Vertex
{
    float x, y, z, w;
    float r, g, b, a;
};
vector< Vertex > verts;
vector< Vertex > verts_bunny;
vector<std::vector<int>> bunny_g;

#define sqr(x) ((x) * (x))

float get_dist(Vertex a, Vertex b) {
    return sqrt(sqr(a.x - b.x) + sqr(a.y - b.y) + sqr(a.z - b.z));
}

struct Triangle {
    Vertex a, b, c;

    Triangle(Vertex a, Vertex b, Vertex c): a(a), b(b), c(c) {}
};
std::vector<Triangle> triangles_bunny;

std::string os = "macOS";
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

void build_triangles_bunny() {
    bunny_g.assign(verts_bunny.size(), std::vector<int>());
    for (int i = 0; i < verts_bunny.size(); ++i) {
        std::cerr << i << std::endl;
        vector<pair<float, int>> ps;
        for (int j = 0; j < verts_bunny.size(); ++j) {
            if (i == j) {
                continue;
            }
            ps.emplace_back(get_dist(verts_bunny[i], verts_bunny[j]), j);
        }
        std::sort(ps.begin(), ps.end());
        for (int j = 0; j < std::min(K, static_cast<int>(ps.size())); ++j) {
            bunny_g[i].push_back(ps[j].second);
        }
    }

    for (int i = 0; i < bunny_g.size(); ++i) {
        for (int j = 0; j < bunny_g[i].size(); ++j) {
            int n1 = bunny_g[i][j];
            for (int k = 0; k < bunny_g[n1].size(); ++k) {
                int n2 = bunny_g[n1][k];
                if (n2 == i || n2 == n1) {
                    continue;
                }
                if (std::find(bunny_g[i].begin(), bunny_g[i].end(), n2) != bunny_g[i].end()) {
                    triangles_bunny.emplace_back(verts_bunny[i], verts_bunny[n1], verts_bunny[n2]);
                }
            }
        }
    }
}

void readVerticesFromPly() {
    std::mt19937 gen;
    std::mt19937::result_type seed = time(0);
    auto real_rand = std::bind(std::uniform_real_distribution<float>(0,1),
                               std::mt19937(seed));

    auto in = std::ifstream("../bunny/bun000.ply");
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
        verts_bunny.push_back(cur);
    }

    build_triangles_bunny();
}

template<class T>
void relax_min(T& a, T b) {
    a = std::min(a, b);
}

template<class T>
void relax_max(T &a, T b) {
    a = std::max(a, b);
}

void print_min_max_coord(std::vector<Vertex> v) {
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

    for (auto tr : triangles_bunny) {
        glBegin(GL_TRIANGLES);

        glColor3f( tr.a.r, tr.a.g, tr.a.b );
        glVertex3f( tr.a.x, tr.a.y, tr.a.z );

        glColor3f( tr.b.r, tr.b.g, tr.b.b );
        glVertex3f( tr.b.x, tr.b.y, tr.b.z );

        glColor3f( tr.c.r, tr.c.g, tr.c.b );
        glVertex3f( tr.c.x, tr.c.y, tr.c.z );

        glEnd();
    }

    // TODO: find nearest points and connect it with each other

    print_min_max_coord(verts);
    print_min_max_coord(verts_bunny);

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
