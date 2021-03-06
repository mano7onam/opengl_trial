#include <GL/glut.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <functional>

struct Vertex
{
    float x, y, z, w;
    float r, g, b, a;
};
std::vector< Vertex > verts;
std::vector< Vertex > verts_bunny;

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

void readVerticesFromPly() {
    std::mt19937 gen;
    std::mt19937::result_type seed = time(0);
    auto real_rand = std::bind(std::uniform_real_distribution<float>(0,1),
                               std::mt19937(seed));

    auto in = std::ifstream("/home/mano/CLionProjects/OpenGL/bunny/bun000.ply");
    std::string dummy;
    for (int i = 0; i < 24; ++i) {
        getline(in, dummy);
    }
    float x, y, z;
    for (int i = 0; i < 40256; ++i) {
        in >> x >> y >> z;
        Vertex cur;
        cur.x = x * 200;
        cur.y = z * 200;
        cur.z = y * 200;
        cur.a = 0.5f;
        cur.r = real_rand();
        cur.g = real_rand();
        cur.b = real_rand();
        verts_bunny.push_back(cur);
    }
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
    gluPerspective( 60.0, w / h, 1.0, 10000.0 );
//    gluPerspective( 60.0, w / h, 0.001, 1.0 );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
//    gluLookAt( 70, 70, 70, 0, 0, 0, 0, 0, 1 );
    gluLookAt( 40, 40, 40, 0, 0, 20, 0, 0, 1 );

    glRotatef( angle, 0, 0, 1 );

    // draw curve
    glEnableClientState( GL_VERTEX_ARRAY );
    glEnableClientState( GL_COLOR_ARRAY );
    glVertexPointer( 3, GL_FLOAT, sizeof( Vertex ), &verts_bunny[0].x );
    glColorPointer( 4, GL_FLOAT, sizeof( Vertex ), &verts_bunny[0].r );
    glDrawArrays( GL_LINE_STRIP, 0, verts_bunny.size() );
    glDisableClientState( GL_VERTEX_ARRAY );
    glDisableClientState( GL_COLOR_ARRAY );

//    for (int i = 0; i < 40000; i += 3) {
//        glBegin(GL_TRIANGLES);
//        for (int j = 0; j < 3; ++j) {
//            auto v = verts[i + j];
//            glColor3f( v.r, v.g, v.b );
//            glVertex3f( v.x, v.y, v.z );
//        }
//        glEnd();
//    }

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