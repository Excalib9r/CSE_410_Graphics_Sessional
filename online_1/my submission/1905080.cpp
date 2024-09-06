#include <bits/stdc++.h>
using namespace std;

#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <glut.h>
#endif

int counter = 0;

void drawCircle(){
    double r = 5;
    double theta = 0;
    double prev_x = 5;
    double prev_y = 0;
    glColor3f(1.0, 1.0, 1.0);
    while(theta  < 360){
        double x = r * cos(theta);
        double y = r * sin(theta);
        glBegin(GL_LINES);
            glVertex2f(prev_x, prev_y);
            glVertex2f(x,y);
        glEnd();
        theta += 1;
        prev_x = x;
        prev_y = y;
    }
}

int animate;

void display()
{
    // glEnable(GL_DEPTH_TEST);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(
        0, 0, 10, 
        0, 0, 0, 
        0, 1, 0
        );
        
    glColor3f(1, 1, 1);
    // square(.25);

    drawCircle();

    // glFlush();
    glutSwapBuffers();
}

void init()
{
    // glClearColor(0.1f, .0f, 0.0f, 1.0f); // Set background color to black and opaque
    animate = 1;


    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, 1, 1, 100);

}


void Timer(int value){
    // printf("We are in Timer function. couter : %d\n", ++counter);
    if(animate){
        counter++;
    }

    glutPostRedisplay();
    glutTimerFunc(10, Timer, 0);
}

int main(int argc, char **argv)
{
    printf("Hello World\n");
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(480, 480);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("OpenGL Demo");

    init();

    glutDisplayFunc(display);

    // glutIdleFunc(idle);
    // glutTimerFunc(0, Timer, 0);
    glutMainLoop();
    return 0;
}