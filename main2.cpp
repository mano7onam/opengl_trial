#include <GL/glut.h>
#include <vector>

#include "ply_loader.h"

float angle = 0;
void timer( int extra )
{
//    // spin
//    angle += 0.5;
//
//    glutPostRedisplay();
//    glutTimerFunc( 16, timer, 0 );
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

    // draw
    Model_PLY *model = new Model_PLY();
    model->Load("/home/mano/CLionProjects/OpenGL/bunny/bun000.ply");
    model->Draw();

    glutSwapBuffers();
}

int main( int argc, char **argv )
{
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );
    glutInitWindowSize( 800, 600 );
    glutCreateWindow( "Bunny" );

    glutDisplayFunc( display );
    glutTimerFunc( 0, timer, 0 );

    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    glEnable( GL_POINT_SMOOTH );
    glPointSize(1.0f);

    glutMainLoop();
}