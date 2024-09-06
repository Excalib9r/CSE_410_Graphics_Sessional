#include <bits/stdc++.h>
#include <windows.h>
#include "bitmap_image.hpp"
#include "1905080.hpp"
#include <GL/glut.h>
using namespace std;

#define pi (2 * acos(0.0))
std::vector<PointLight> pointLights;
std::vector<SpotLight> spotLights;
class Material;
std::vector<Material *> materials;
int recursionLevel;

Point up(0, 0, 1);
Point rightSide(-1 / sqrt(2), 1 / sqrt(2), 0);
Point direction(-1 / sqrt(2), -1 / sqrt(2), 0);

Point eye(150, 0, 20);

bitmap_image image;
int imageCount = 1;

double angleZ = pi / 180;
int numSegments;

int imageHeight, imageWidth;

double windowWidth = 500, windowHeight = 500;
double viewAngle = 80;

void setColorValues(int nearestObjectIndex, Raylight ray, double &color_x, double &color_y, double &color_z)
{
    color_x = 0;
    color_y = 0;
    color_z = 0;
    double t = materials[nearestObjectIndex]->intersect(ray, color_x, color_y, color_z, 1);
    color_x = min(max(color_x, 0.0), 1.0);
    color_y = min(max(color_y, 0.0), 1.0);
    color_z = min(max(color_z, 0.0), 1.0);
}

int findNearestObjectIndex(Raylight ray, double &tMin, double &color_x, double &color_y, double &color_z)
{
    int nearestObjectIndex = -1;
    double t;
    tMin = -1;
    for (int k = 0; k < (int)materials.size(); k++)
    {
        t = materials[k]->intersect(ray, color_x, color_y, color_z, 0);
        if (t > 0 && (nearestObjectIndex == -1 || t < tMin))
            tMin = t, nearestObjectIndex = k;
    }
    return nearestObjectIndex;
}

double calculatePlaneDistance()
{
    return (windowHeight / 2.0) / tan((pi * viewAngle) / (360.0));
}

Point calculateTopLeft(double planeDistance)
{
    Point topLeft = addPoint(eye, addPoint(mulConstant(direction, planeDistance), subPoint(mulConstant(up, (windowHeight / 2.0)), mulConstant(rightSide, (windowWidth / 2.0)))));
    double du = windowWidth * 1.0 / (imageWidth);
    double dv = windowHeight * 1.0 / (imageHeight);
    return addPoint(topLeft, subPoint(mulConstant(rightSide, du / 2.0), mulConstant(up, dv / 2.0)));
}

void capture()
{
    cout << "Image Captured" << endl;

    for (int i = 0; i < imageWidth; i++)
        for (int j = 0; j < imageHeight; j++)
            image.set_pixel(i, j, 0, 0, 0);

    double planeDistance = calculatePlaneDistance();
    Point topLeft = calculateTopLeft(planeDistance);
    int nearestObjectIndex = -1;
    double t, tMin;

    double du = windowWidth * 1.0 / (imageWidth);
    double dv = windowHeight * 1.0 / (imageHeight);

    for (int i = 0; i < imageWidth; i++)
    {
        for (int j = 0; j < imageHeight; j++)
        {
            Point pixel = addPoint(topLeft, subPoint(mulConstant(rightSide, du * i), mulConstant(up, dv * j)));
            Raylight ray(eye, subPoint(pixel, eye));
            double color_x;
            double color_y;
            double color_z;

            nearestObjectIndex = findNearestObjectIndex(ray, tMin, color_x, color_y, color_z);

            if (nearestObjectIndex != -1)
            {
                setColorValues(nearestObjectIndex, ray, color_x, color_y, color_z);

                image.set_pixel(i, j, 255 * color_x, 255 * color_y, 255 * color_z);
            }
        }
    }

    image.save_image("Image_" + to_string(imageCount) + ".bmp");
    imageCount++;
    cout << "Image saved" << endl;
}

void axes()
{
    glLineWidth(1);
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

void grid()
{
    glColor3f(0.6, 0.6, 0.6);
    glBegin(GL_LINES);
    {
        for (int i = -8; i <= 8; i++)
        {
            if (i == 0)
                continue;
            glVertex3f(i * 10, -90, 0);
            glVertex3f(i * 10, 90, 0);

            glVertex3f(-90, i * 10, 0);
            glVertex3f(90, i * 10, 0);
        }
    }
    glEnd();
}

Material *createSphere(ifstream &in)
{
    double x, y, z;
    in >> x >> y >> z;
    Point ref_point(x, y, z);
    double length;
    in >> length;
    double color_x, color_y, color_z;
    in >> color_x >> color_y >> color_z;
    double ambient, diffuse, specular, reflection;
    in >> ambient >> diffuse >> specular >> reflection;
    int shininess;
    in >> shininess;
    Sphere *obj = new Sphere(ref_point, length);
    obj->setColor(color_x, color_y, color_z);
    obj->setAmbient(ambient);
    obj->setDiffuse(diffuse);
    obj->setSpecular(specular);
    obj->setReflection(reflection);
    obj->setShininess(shininess);
    return obj;
}

Material *createTriangle(ifstream &in)
{
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    in >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
    Point A(x1, y1, z1);
    Point B(x2, y2, z2);
    Point C(x3, y3, z3);
    double color_x, color_y, color_z;
    in >> color_x >> color_y >> color_z;
    double ambient, diffuse, specular, reflection;
    in >> ambient >> diffuse >> specular >> reflection;
    int shininess;
    in >> shininess;
    Triangle *obj = new Triangle(A, B, C);
    obj->setColor(color_x, color_y, color_z);
    obj->setAmbient(ambient);
    obj->setDiffuse(diffuse);
    obj->setSpecular(specular);
    obj->setReflection(reflection);
    obj->setShininess(shininess);
    return obj;
}

Material *createGeneralObject(ifstream &in)
{
    double A, B, C, D, E, F, G, H, I, J;
    in >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
    double x, y, z;
    in >> x >> y >> z;
    Point ref_point(x, y, z);
    double height, width, length;
    in >> height >> width >> length;
    double color_x, color_y, color_z;
    in >> color_x >> color_y >> color_z;
    double ambient, diffuse, specular, reflection;
    in >> ambient >> diffuse >> specular >> reflection;
    int shininess;
    in >> shininess;
    GenaeralObject *obj = new GenaeralObject(A, B, C, D, E, F, G, H, I, J);
    obj->setColor(color_x, color_y, color_z);
    obj->setAmbient(ambient);
    obj->setDiffuse(diffuse);
    obj->setSpecular(specular);
    obj->setReflection(reflection);
    obj->setShininess(shininess);
    obj->ref_point = ref_point;
    return obj;
}

void loadData()
{
    ifstream in("scene.txt");
    in >> recursionLevel >> imageHeight;

    imageWidth = imageHeight;

    int objCount;
    in >> objCount;

    for (int i = 0; i < objCount; i++)
    {
        string objType;
        in >> objType;

        Material *obj;

        if (objType == "sphere")
        {
            obj = createSphere(in);
        }
        else if (objType == "triangle")
        {
            obj = createTriangle(in);
        }
        else if (objType == "general")
        {
            obj = createGeneralObject(in);
        }
        else
        {
            cout << objType << " is not a valid object type" << endl;
            continue;
        }
        materials.push_back(obj);
    }

    int lightCount;
    in >> lightCount;

    for (int i = 0; i < lightCount; i++)
    {
        double x, y, z;
        in >> x >> y >> z;
        Point position(x, y, z);
        double color_x, color_y, color_z;
        in >> color_x >> color_y >> color_z;
        PointLight light = PointLight(position, color_x, color_y, color_z);
        pointLights.push_back(light);
    }

    int spotlightCount;
    in >> spotlightCount;

    for (int i = 0; i < spotlightCount; i++)
    {
        double x, y, z;
        in >> x >> y >> z;
        Point position = Point(x, y, z);
        double color_x, color_y, color_z;
        in >> color_x >> color_y >> color_z;
        double dir_x, dir_y, dir_z;
        in >> dir_x >> dir_y >> dir_z;
        Point direction = Point(dir_x, dir_y, dir_z);
        double cutoffAngle;
        in >> cutoffAngle;
        SpotLight spotlight = SpotLight(position, direction, color_x, color_y, color_z, cutoffAngle);
        spotLights.push_back(spotlight);
    }

    Material *floor;
    floor = new Floor(400, 10);
    floor->setColor(0.5, 0.5, 0.5);
    double ambient = 0.4;
    double diffuse = 0.2;
    double specular = 0.2;
    double reflection = 0.2;
    floor->setAmbient(ambient);
    floor->setDiffuse(diffuse);
    floor->setSpecular(specular);
    floor->setReflection(reflection);
    materials.push_back(floor);
}

void init()
{
    numSegments = 36;

    loadData();
    image = bitmap_image(imageWidth, imageHeight);

    glClearColor(0, 0, 0, 0);

    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();

    gluPerspective(80, 1, 1, 1000.0);
}
void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye.x, eye.y, eye.z,
              eye.x + direction.x, eye.y + direction.y, eye.z + direction.z,
              up.x, up.y, up.z);
    glMatrixMode(GL_MODELVIEW);

    axes();
    grid();

    for (int i = 0; i < materials.size(); i++)
    {
        Material *object = materials[i];
        object->draw();
    }

    for (int i = 0; i < pointLights.size(); i++)
    {
        pointLights[i].draw();
    }

    for (int i = 0; i < spotLights.size(); i++)
    {
        spotLights[i].draw();
    }

    glutSwapBuffers();
}

void animate()
{
    // codes for any changes in Models, Camera
    glutPostRedisplay();
}

void keyboardHandler(unsigned char key, int x, int y)
{
    double view = 0.2;
    double rate = angleZ;
    switch (key)
    {
    case 'w':
        eye.z += view;
        break;
    case 's':
        eye.z -= view;
        break;
    case 'd':
        angleZ -= 3.0f;
        break;
    case 'a':
        angleZ += 3.0f;
        break;
    case '0':
        capture();
        break;
    case '1':
        rightSide = addPoint(mulConstant(rightSide, cos(rate)), mulConstant(direction, sin(rate)));
        direction = subPoint(mulConstant(direction, cos(rate)), mulConstant(rightSide, sin(rate)));
        break;
    case '2':
        rightSide = addPoint(mulConstant(rightSide, cos(-rate)), mulConstant(direction, sin(-rate)));
        direction = subPoint(mulConstant(direction, cos(-rate)), mulConstant(rightSide, sin(-rate)));
        break;
    case '3':
        direction = addPoint(mulConstant(direction, cos(rate)), mulConstant(up, sin(rate)));
        up = subPoint(mulConstant(up, cos(rate)), mulConstant(direction, sin(rate)));
        break;
    case '4':
        direction = addPoint(mulConstant(direction, cos(-rate)), mulConstant(up, sin(-rate)));
        up = subPoint(mulConstant(up, cos(-rate)), mulConstant(direction, sin(-rate)));
        break;
    case '5':
        up = addPoint(mulConstant(up, cos(rate)), mulConstant(rightSide, sin(rate)));
        rightSide = subPoint(mulConstant(rightSide, cos(rate)), mulConstant(up, sin(rate)));
        break;
    case '6':
        up = addPoint(mulConstant(up, cos(-rate)), mulConstant(rightSide, sin(-rate)));
        rightSide = subPoint(mulConstant(rightSide, cos(-rate)), mulConstant(up, sin(-rate)));
        break;
    default:
        break;
    }
    glutPostRedisplay();
}
void specialKey(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_UP:
        eye.x = eye.x + direction.x;
        eye.y = eye.y + direction.y;
        eye.z = eye.z + direction.z;
        ;
        break;
    case GLUT_KEY_DOWN:
        eye.x = eye.x - direction.x;
        eye.y = eye.y - direction.y;
        eye.z = eye.z - direction.z;
        break;
    case GLUT_KEY_RIGHT:
        eye.x = eye.x + rightSide.x;
        eye.y = eye.y + rightSide.y;
        eye.z = eye.z + rightSide.z;
        break;
    case GLUT_KEY_LEFT:
        eye.x = eye.x - rightSide.x;
        eye.y = eye.y - rightSide.y;
        eye.z = eye.z - rightSide.z;
        break;
    case GLUT_KEY_PAGE_UP:
        eye.x = eye.x + up.x;
        eye.y = eye.y + up.y;
        eye.z = eye.z + up.z;
        break;
    case GLUT_KEY_PAGE_DOWN:
        eye.x = eye.x - up.x;
        eye.y = eye.y - up.y;
        eye.z = eye.z - up.z;
        break;
    default:
        break;
    }
    glutPostRedisplay();
}
int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

    glutCreateWindow("Ray Tracing");

    init();

    glEnable(GL_DEPTH_TEST);

    glutDisplayFunc(display);
    glutIdleFunc(animate);

    glutKeyboardFunc(keyboardHandler);
    glutSpecialFunc(specialKey);

    glutMainLoop();
    return 0;
}