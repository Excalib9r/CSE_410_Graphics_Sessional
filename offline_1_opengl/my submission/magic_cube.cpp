#include <bits/stdc++.h>

#ifdef __linux__
    #include<GL/glut.h>
#elif WIN32
    #include <windows.h>
    #include <GL/glut.h>
#endif

class Point
{
public:
    double x, y, z;
    Point() {
        x = 0;
        y = 0;
        z = 0;
    }
    Point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    double getLength(){
        return sqrt(x*x+y*y+z*z);
    }
};

Point normalPoint(Point a){
    double length = a.getLength();
    return Point(a.x/length, a.y/length, a.z/length);
}

Point mulConstant(Point a, double c){
    return Point(a.x*c, a.y*c, a.z*c);
}

Point addPoint(Point a, Point b){
    return Point(a.x+b.x, a.y+b.y, a.z+b.z);
}

Point subPoint(Point a, Point b){
    return Point(a.x-b.x, a.y-b.y, a.z-b.z);
}

class Sphere{
    public:
    int numPoints;
    double radius;
    Point **points;
    Sphere(int numPoints, double radius){
        this->numPoints = numPoints;
        this->radius = radius;
        points = new Point*[numPoints];
        for(int i = 0; i < numPoints; i++){
            points[i] = new Point[numPoints];
        }
    }
    void setNumPoints(int numPoints){
        this->numPoints = numPoints;
    }
    void setRadius(double radius){
        this->radius = radius;
    }
    void generatePoints(){
        double x,y;
        for(int i = 0;i < numPoints; i++){
            for(int j = 0; j < numPoints; j++){
                x = -1 + (double)i/numPoints*2;
                y = -1 + (double)j/numPoints*2;
                points[i][j].x=x;
                points[i][j].y=y;
                points[i][j].z=1;

                points[i][j] = normalPoint(points[i][j]);
                points[i][j] = mulConstant(points[i][j], radius);
            }
        }
    }
    void draw(){
        for(int i = 0; i< numPoints - 1; i++){
            for(int j = 0; j < numPoints - 1; j++){
                glBegin(GL_QUADS);{
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                }glEnd();
            }
	    }
    }
};

class Cylinder{
    double depth;
    double offset = 70.5287794*M_PI/180.0;
    double radius;
    int numPoints;
    Point *points;
    public:
    Cylinder(double depth, double radius, int numPoints){
        this->depth = depth;
        this->radius = radius;
        this->numPoints = numPoints;
        points = new Point[numPoints];
    }
    void setDepth(double depth){
        this->depth = depth;
    }
    void setRadius(double radius){
        this->radius = radius;
    }
    void setNumPoints(int numPoints){
        this->numPoints = numPoints;
    }
    void generatePoints(){
        for (int i = 0; i < numPoints; i++) {
            double angle = -offset/2 +  i * offset / numPoints;
            points[i].x = radius * cos(angle);
            points[i].y = radius * sin(angle);
        }
    }
    void draw(){
        glBegin(GL_QUADS);
            for (int i = 0; i < numPoints - 1; i++) {
                glVertex3f(points[i].x, points[i].y, depth/2);
                glVertex3f(points[i].x, points[i].y, -depth/2);
                glVertex3f(points[i+1].x, points[i+1].y, -depth/2);
                glVertex3f(points[i+1].x, points[i+1].y, depth/2);
            }
        glEnd();
    }
};

Point up(0, 0, 1);
Point rightSide(0, 1, 0);
Point direction(-1, 0, 0);

double eye_x = 6, eye_y = 0, eye_z = 1;

double maxLength = 2.5;
double tLength = 2.5;
double maxRadious = maxLength / sqrt(3.0);
double sRadious = 0;
double step = maxRadious/25.0;

double angleZ = 30.0;

Sphere sphere(101, sRadious);
Cylinder cylinder(tLength*sqrt(2), sRadious, 101);

void axes() {
    glLineWidth(2);
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
        glVertex3f(100, 0, 0);
        glVertex3f(-100, 0, 0);

        glVertex3f(0, -100, 0);
        glVertex3f(0, 100, 0);

        glVertex3f(0, 0, 100);
        glVertex3f(0, 0, -100);
    glEnd();
}
void triangle()
{
	glBegin(GL_TRIANGLES);{
		glVertex3f(1,0,0);
		glVertex3f(0,1,0);
		glVertex3f(0,0,1);
	}glEnd();
}
void cylinderHelper(){
    glTranslatef(tLength/sqrt(2),0,0);
    cylinder.setDepth(tLength*sqrt(2));
    cylinder.setRadius(sRadious);
    cylinder.generatePoints();
    cylinder.draw();
}
void cylinders(){
    glColor3f(0.8f, 1.0f, 0.0f);
    for(int i=0;i<4;i++){
        glPushMatrix();{
            glRotatef(45+i*90,0,1,0);
            cylinderHelper();
        }glPopMatrix();     
    }
    for(int i=0;i<4;i++){
        glPushMatrix();{
            glRotatef(90,1,0,0);
            glRotatef(45+i*90,0,1,0);
            cylinderHelper();
        }glPopMatrix(); 
    }
    for(int i=0;i<4;i++){
        glPushMatrix();{
            glRotatef(90,0,0,1);
            glRotatef(45+i*90,0,1,0);
            cylinderHelper();
        }glPopMatrix();
    }
}
void sphereHelper(int i, int j){
    glPushMatrix();{
        if(j == 0){
             glColor3f(0.3f, (i%2)*0.3, (i+1)%2);
            glRotatef(90*i,0,1,0);
        }
        else{
            glColor3f(0.5f, 0.2f, 0.0f); 
            glRotatef(90+180*i,1,0,0);
        }
        glTranslatef(0,0,tLength);
        sphere.setRadius(sRadious);
        sphere.generatePoints();
        sphere.draw();
        }glPopMatrix();
}
void spheres(){
    for(int i=0;i<4;i++){    
        sphereHelper(i, 0);
    }
    for(int i=0;i<2;i++){    
        sphereHelper(i, 1);
    }
}
void pyramidsHelper(double scale, int i, int j){
    glPushMatrix();{
        if(j == 0){
            glColor3f((i+1)%2, (i%2)*0.5, 0.0f);
            glRotatef(90*i,0,1,0);
        }
        else{
            glColor3f((i%2)*0.5, (i+1)%2, 0.0f);
            glRotatef(90*i,0,1,0);
            glRotatef(180,1,0,1);
        }
        glTranslatef(scale,scale,scale);
        glScaled(tLength,tLength,tLength);
        triangle();
    }glPopMatrix();
}
void pyramids(){
    double scale = maxLength - tLength;
    scale = scale/3.0;
    for(int i=0;i<4;i++){
        pyramidsHelper(scale, i, 0);
    }
    for(int i=0;i<4;i++){
        pyramidsHelper(scale, i, 1);
    }
}
void init()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, 1, 1, 100);
}
void display() {
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);       
    glLoadIdentity();                      
   
    gluLookAt(eye_x,eye_y,eye_z, 
    eye_x+direction.x,eye_y+direction.y,eye_z+direction.z,
    up.x,up.y,up.z);

    glRotatef(angleZ, 0,0 , 1);
    axes();
    pyramids();
    spheres();
    cylinders();

    glutSwapBuffers();
}
void keyboardHandler(unsigned char key, int x, int y) {
    double view = 0.2;
    double rate = 0.1;
    switch (key) {
        case ',':
            tLength-=0.1;
            sRadious+=step;
            if( tLength < 0 ) {
                tLength = 0;
                sRadious = maxRadious;
            }
            break;
        case '.':
            tLength+=0.1;
            sRadious-=step;
            if( tLength > maxLength ) {
                tLength = maxLength;
                sRadious = 0;
            }
            break;
        case 'w':
            eye_z += view;
            break;
        case 's':
            eye_z -= view;
            break;
        case 'd':
            angleZ -= 3.0f;
            break;
        case 'a':
            angleZ += 3.0f;
            break;
        case '1':
            rightSide = addPoint(mulConstant(rightSide,cos(rate)), mulConstant(direction,sin(rate)));
            direction = subPoint(mulConstant(direction,cos(rate)),mulConstant(rightSide,sin(rate)));
            break;
        case '2':
            rightSide = addPoint(mulConstant(rightSide,cos(-rate)), mulConstant(direction,sin(-rate)));
            direction = subPoint(mulConstant(direction,cos(-rate)),mulConstant(rightSide,sin(-rate)));
            break;
        case '3':
            direction = addPoint(mulConstant(direction,cos(rate)), mulConstant(up,sin(rate)));
            up = subPoint(mulConstant(up,cos(rate)),mulConstant(direction,sin(rate)));
            break;
        case '4':
            direction = addPoint(mulConstant(direction,cos(-rate)), mulConstant(up,sin(-rate)));
            up = subPoint(mulConstant(up,cos(-rate)),mulConstant(direction,sin(-rate)));
            break;
        case '5':
            up = addPoint(mulConstant(up,cos(rate)),mulConstant(rightSide,sin(rate)));
            rightSide = subPoint(mulConstant(rightSide,cos(rate)), mulConstant(up,sin(rate)));
            break;
        case '6':
            up = addPoint(mulConstant(up,cos(-rate)),mulConstant(rightSide,sin(-rate)));
            rightSide = subPoint(mulConstant(rightSide,cos(-rate)), mulConstant(up,sin(-rate)));
            break;
       default:
            break;
    }
    glutPostRedisplay();
}
void specialKey(int key, int x,int y)
{
	switch(key){
		case GLUT_KEY_UP:		
			eye_x = eye_x + direction.x;
            eye_y = eye_y + direction.y;
            eye_z = eye_z + direction.z;;
			break;
		case GLUT_KEY_DOWN:		
			eye_x = eye_x - direction.x;
            eye_y = eye_y - direction.y;
            eye_z = eye_z - direction.z;
			break;
		case GLUT_KEY_RIGHT:
			eye_x = eye_x + rightSide.x;
            eye_y = eye_y + rightSide.y;
            eye_z = eye_z + rightSide.z;
			break;
		case GLUT_KEY_LEFT :
			eye_x = eye_x - rightSide.x;
            eye_y = eye_y - rightSide.y;
            eye_z = eye_z - rightSide.z;
			break;
		case GLUT_KEY_PAGE_UP:
            eye_x = eye_x + up.x;
            eye_y = eye_y + up.y;
            eye_z = eye_z + up.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
            eye_x = eye_x - up.x;
            eye_y = eye_y - up.y;
            eye_z = eye_z - up.z;
			break;
		default:
			break;
	}
	glutPostRedisplay();
}
int main(int argc,char** argv){
    glutInit(&argc,argv);
    glutInitWindowSize(450, 450);
    glutInitWindowPosition(750, 250);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Test");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboardHandler);
    glutSpecialFunc(specialKey);
    init();
    glutMainLoop();
    return 0;
}