#include <GL/glut.h>
#include <vector>

struct Vertex
{
    float x, y, z, w;
    float r, g, b, a;
};
std::vector< Vertex > verts;

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

void display(void)
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    const double w = glutGet( GLUT_WINDOW_WIDTH );
    const double h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 60.0, w / h, 1.0, 10000.0 );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt( 70, 70, 70, 0, 0, 0, 0, 0, 1 );

    glRotatef( angle, 0, 0, 1 );

    // draw curve
    glEnableClientState( GL_VERTEX_ARRAY );
    glEnableClientState( GL_COLOR_ARRAY );
    glVertexPointer( 3, GL_FLOAT, sizeof( Vertex ), &verts[0].x );
    glColorPointer( 4, GL_FLOAT, sizeof( Vertex ), &verts[0].r );
    glDrawArrays( GL_LINE_STRIP, 0, verts.size() );
    glDisableClientState( GL_VERTEX_ARRAY );
    glDisableClientState( GL_COLOR_ARRAY );

    glutSwapBuffers();
}

int main( int argc, char **argv )
{
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );
    glutInitWindowSize( 800,600 );
    glutCreateWindow( "Attractor" );

    glutDisplayFunc( display );
    glutTimerFunc( 0, timer, 0 );

    fillVerts();

    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    glEnable( GL_POINT_SMOOTH );
    glPointSize(1.0f);

    glutMainLoop();
    return 0;
}