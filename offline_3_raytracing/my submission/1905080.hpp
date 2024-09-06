#include <bits/stdc++.h>
#include "bitmap_image.hpp"
#include <GL/glut.h>

using namespace std;

#define pi (2 * acos(0.0))
extern bitmap_image image;

class Point
{
public:
    double x, y, z;
    Point()
    {
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
    double getLength()
    {
        return sqrt(x * x + y * y + z * z);
    }
};

double mulPoint(Point a, Point b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Point normalPoint(Point a)
{
    double length = a.getLength();
    return Point(a.x / length, a.y / length, a.z / length);
}

Point mulConstant(Point a, double c)
{
    return Point(a.x * c, a.y * c, a.z * c);
}

Point addPoint(Point a, Point b)
{
    return Point(a.x + b.x, a.y + b.y, a.z + b.z);
}

Point subPoint(Point a, Point b)
{
    return Point(a.x - b.x, a.y - b.y, a.z - b.z);
}

Point crossProduct(Point a, Point b)
{
    return Point(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

class PointLight
{
public:
    Point position;
    double color_x, color_y, color_z;
    PointLight(Point position, double color_x, double color_y, double color_z)
    {
        this->position = position;
        this->color_x = color_x;
        this->color_y = color_y;
        this->color_z = color_z;
    }
    void draw()
    {
        glPointSize(5);
        glBegin(GL_POINTS);
        glColor3f(color_x, color_y, color_z);
        glVertex3f(position.x, position.y, position.z);
        glEnd();
    }
};

class SpotLight
{
public:
    Point position;
    Point direction;
    double color_x, color_y, color_z;
    double cutoffAngle;
    SpotLight(Point position, Point direction, double color_x, double color_y, double color_z, double cutoffAngle)
    {
        this->position = position;
        this->direction = direction;
        this->color_x = color_x;
        this->color_y = color_y;
        this->color_z = color_z;
        this->cutoffAngle = cutoffAngle;
    }
    void draw()
    {
        glPointSize(15);
        glBegin(GL_POINTS);
        glColor3f(color_x, color_y, color_z);
        glVertex3f(position.x, position.y, position.z);
        glEnd();
    }
};

class Raylight
{
public:
    Point position;
    Point normal;
    Raylight()
    {
        this->position = Point(0, 0, 0);
        this->normal = Point(0, 0, 0);
    }
    Raylight(Point position, Point direction)
    {
        this->position = position;
        this->normal = normalPoint(direction);
    }
};

class Material;

extern std::vector<PointLight> pointLights;
extern std::vector<SpotLight> spotLights;
extern std::vector<Material *> materials;
extern int recursionLevel;

void calculateColorComponent(double &color, double lightColor, double diffuse, double val, double colorAtIntersection, double specular, double phong, double shininess)
{
    color += lightColor * diffuse * val * colorAtIntersection;
    color += lightColor * specular * pow(phong, shininess) * colorAtIntersection;
}

double calculateColorComponent2(double lightColor, double diffuse, double val, double colorAtIntersection, double specular, double phong, double shininess)
{
    double color = lightColor * diffuse * val * colorAtIntersection;
    color += lightColor * specular * pow(phong, shininess) * colorAtIntersection;
    return color;
}

class Material
{
public:
    double ambient, diffuse, specular, reflection;
    int shininess;
    Point ref_point;
    double height, width, length;
    double color_x, color_y, color_z;
    Material()
    {
        this->ambient = 0;
        this->diffuse = 0;
        this->specular = 0;
        this->reflection = 0;
    }
    Material(Point ref_point, double height, double width, double length, double color_x, double color_y, double color_z, double ambient, double diffuse, double specular, double reflection, int shininess)
    {
        this->ref_point = ref_point;
        this->height = height;
        this->width = width;
        this->length = length;
        this->color_x = color_x;
        this->color_y = color_y;
        this->color_z = color_z;
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->reflection = reflection;
        this->shininess = shininess;
    }
    void setColor(double color_x, double color_y, double color_z)
    {
        this->color_x = color_x;
        this->color_y = color_y;
        this->color_z = color_z;
    }
    virtual void getColor(double &color_x, double &color_y, double &color_z, Point point)
    {
        color_x = this->color_x;
        color_y = this->color_y;
        color_z = this->color_z;
    }

    void setShininess(int shininess)
    {
        this->shininess = shininess;
    }
    void setAmbient(double ambient)
    {
        this->ambient = ambient;
    }
    void setDiffuse(double diffuse)
    {
        this->diffuse = diffuse;
    }
    void setSpecular(double specular)
    {
        this->specular = specular;
    }
    void setReflection(double reflection)
    {
        this->reflection = reflection;
    }
    virtual void draw() = 0;
    virtual double helperFunc(Raylight ray, double &color_x, double &color_y, double &color_z, int level) = 0;
    virtual Raylight getNormal(Point point, Raylight incident) = 0;
    bool getObsecured(Point intersection, int i, Raylight &lightRay, Raylight &normalRay, double &color_x, double &color_y, double &color_z)
    {
        bool obscured = false;
        Point positionOfLight = pointLights[i].position;
        Point directionOfLight = subPoint(intersection, positionOfLight);
        Point normal = normalPoint(directionOfLight);

        lightRay = Raylight(positionOfLight, normal);
        normalRay = getNormal(intersection, lightRay);

        double t2 = subPoint(intersection, positionOfLight).getLength();
        if (t2 < 1e-5)
            return true;

        for (Material *material : materials)
        {
            double t3 = material->helperFunc(lightRay, color_x, color_y, color_z, 0);
            if (t3 > 0 && t3 + 1e-5 < t2)
            {
                obscured = true;
                break;
            }
        }
        return obscured;
    }

    Raylight calculateReflection(Point intersection, Raylight lightRay, Raylight normalRay)
    {
        double helper2 = mulPoint(lightRay.normal, normalRay.normal);
        Point reflectionHelper = mulConstant(normalRay.normal, 2 * helper2);
        return Raylight(intersection, subPoint(lightRay.normal, reflectionHelper));
    }

    bool isObscured(double t2, Raylight lightRay, double color_x, double color_y, double color_z)
    {
        for (Material *obj : materials)
        {
            double t3 = obj->helperFunc(lightRay, color_x, color_y, color_z, 0);
            if (t3 > 0 && t3 + 1e-5 < t2)
            {
                return true;
            }
        }
        return false;
    }

    virtual double intersect(Raylight ray, double &color_x, double &color_y, double &color_z, int level)
    {
        double t = helperFunc(ray, color_x, color_y, color_z, level);

        if (t < 0)
            return -1;
        if (level == 0)
            return t;

        Point intersection = addPoint(ray.position, mulConstant(ray.normal, t));
        double colorAtIntersection[3];
        getColor(colorAtIntersection[0], colorAtIntersection[1], colorAtIntersection[2], intersection);

        color_x = colorAtIntersection[0] * ambient;
        color_y = colorAtIntersection[1] * ambient;
        color_z = colorAtIntersection[2] * ambient;

        for (int i = 0; i < pointLights.size(); i++)
        {
            Raylight lightRay, normalRay;
            bool obscured = getObsecured(intersection, i, lightRay, normalRay, color_x, color_y, color_z);
            if (!obscured)
            {
                double value = mulPoint(lightRay.normal, normalRay.normal);
                double val = max(0.0, -value);

                Point helper1 = mulConstant(normalRay.normal, 2 * value);
                Point helper2 = subPoint(lightRay.normal, helper1);
                Raylight reflection = Raylight(intersection, helper2);
                double value2 = mulPoint(reflection.normal, ray.normal);
                double phong = max(0.0, -value2);

                calculateColorComponent(color_x, pointLights[i].color_x, diffuse, val, colorAtIntersection[0], specular, phong, shininess);
                calculateColorComponent(color_y, pointLights[i].color_y, diffuse, val, colorAtIntersection[1], specular, phong, shininess);
                calculateColorComponent(color_z, pointLights[i].color_z, diffuse, val, colorAtIntersection[2], specular, phong, shininess);
            }
        }

        for (int i = 0; i < spotLights.size(); i++)
        {

            Point lightPosition = spotLights[i].position;
            Point lightDirection = subPoint(intersection, lightPosition);
            Point normal = normalPoint(lightDirection);

            double dot = mulPoint(lightDirection, spotLights[i].direction);
            double helper1 = lightDirection.getLength() * spotLights[i].direction.getLength();
            double angle = acos(dot / helper1) * (180.0 / pi);

            if (fabs(angle) < spotLights[i].cutoffAngle)
            {
                Raylight lightRay = Raylight(lightPosition, lightDirection);
                Raylight normalRay = getNormal(intersection, lightRay);

                Raylight reflection = calculateReflection(intersection, lightRay, normalRay);

                double t2 = subPoint(intersection, lightPosition).getLength();
                if (t2 < 1e-5)
                    continue;

                bool obscured = false;

                for (Material *obj : materials)
                {
                    double t3 = obj->helperFunc(lightRay, color_x, color_y, color_z, 0);
                    if (t3 > 0 && t3 + 1e-5 < t2)
                    {
                        obscured = true;
                        break;
                    }
                }

                if (!obscured)
                {
                    double value = mulPoint(ray.normal, reflection.normal);
                    double phong = max(0.0, -value);
                    double value2 = mulPoint(lightRay.normal, normalRay.normal);
                    double val = max(0.0, -value2);

                    color_x = calculateColorComponent2(spotLights[i].color_x, diffuse, val, colorAtIntersection[0], specular, phong, shininess);
                    color_y = calculateColorComponent2(spotLights[i].color_y, diffuse, val, colorAtIntersection[1], specular, phong, shininess);
                    color_z = calculateColorComponent2(spotLights[i].color_z, diffuse, val, colorAtIntersection[2], specular, phong, shininess);
                }
            }
        }
        if (level < recursionLevel)
        {
            Raylight normalRay = getNormal(intersection, ray);

            double helper2 = mulPoint(ray.normal, normalRay.normal);
            Point reflectionHelper = mulConstant(normalRay.normal, 2 * helper2);
            Raylight reflectionRay = Raylight(intersection, subPoint(ray.normal, reflectionHelper));

            reflectionRay.position = addPoint(reflectionRay.position, mulConstant(reflectionRay.normal, 1e-5));

            int closestObject = -1;
            double t = -1, tMin = 1e9;

            for (int k = 0; k < (int)materials.size(); k++)
            {
                t = materials[k]->intersect(reflectionRay, color_x, color_y, color_z, 0);
                if (t > 0 && t < tMin)
                    tMin = t, closestObject = k;
            }

            if (closestObject != -1)
            {

                double tempColor_x = 0, tempColor_y = 0, tempColor_z = 0;

                double t = materials[closestObject]->intersect(reflectionRay, tempColor_x, tempColor_y, tempColor_z, level + 1);

                color_x += tempColor_x * reflection;
                color_y += tempColor_y * reflection;
                color_z += tempColor_z * reflection;
            }
        }

        return t;
    }
};

bool checkDimension(double dimensionSize, double pointValue, double refPointValue)
{
    if (fabs(dimensionSize) > 1e-5)
    {
        if (pointValue < refPointValue || pointValue > refPointValue + dimensionSize)
        {
            return false;
        }
    }
    return true;
}

class GenaeralObject : public Material
{
public:
    double A, B, C, D, E, F, G, H, I, J;
    GenaeralObject(double A, double B, double C, double D, double E, double F, double G, double H, double I, double J)
    {
        this->A = A;
        this->B = B;
        this->C = C;
        this->D = D;
        this->E = E;
        this->F = F;
        this->G = G;
        this->H = H;
        this->I = I;
        this->J = J;
    }

    virtual void draw()
    {
        return;
    }
    virtual Raylight getNormal(Point point, Raylight incidentRay)
    {
        Point dir(2 * A * point.x + D * point.y + E * point.z + G,
                  2 * B * point.y + D * point.x + F * point.z + H,
                  2 * C * point.z + E * point.x + F * point.y + I);

        return Raylight(point, dir);
    }

    bool isIntersecting(Point point)
    {
        if (!checkDimension(length, point.x, ref_point.x))
        {
            return false;
        }

        if (!checkDimension(width, point.y, ref_point.y))
        {
            return false;
        }

        if (!checkDimension(height, point.z, ref_point.z))
        {
            return false;
        }

        return true;
    }
    std::tuple<double, double, double> calculateCoefficients(Raylight ray)
    {
        double X0 = ray.position.x, Y0 = ray.position.y, Z0 = ray.position.z;
        double X1 = ray.normal.x, Y1 = ray.normal.y, Z1 = ray.normal.z;

        double C0 = A * X1 * X1 + B * Y1 * Y1 + C * Z1 * Z1 + D * X1 * Y1 + E * X1 * Z1 + F * Y1 * Z1;
        double C1 = 2 * A * X0 * X1 + 2 * B * Y0 * Y1 + 2 * C * Z0 * Z1 + D * (X0 * Y1 + X1 * Y0) + E * (X0 * Z1 + X1 * Z0) + F * (Y0 * Z1 + Y1 * Z0) + G * X1 + H * Y1 + I * Z1;
        double C2 = A * X0 * X0 + B * Y0 * Y0 + C * Z0 * Z0 + D * X0 * Y0 + E * X0 * Z0 + F * Y0 * Z0 + G * X0 + H * Y0 + I * Z0 + J;

        return {C0, C1, C2};
    }

    double calculateDiscriminant(double C0, double C1, double C2)
    {
        return C1 * C1 - 4 * C0 * C2;
    }

    std::pair<double, double> calculateRoots(double C0, double C1, double discriminant)
    {
        double t1 = (-C1 - sqrt(discriminant)) / (2 * C0);
        double t2 = (-C1 + sqrt(discriminant)) / (2 * C0);

        if (t2 < t1)
            std::swap(t1, t2);

        return {t1, t2};
    }

    double findIntersection(Raylight ray, double t1, double t2)
    {
        if (t1 > 0)
        {
            Point intersectionPoint = addPoint(ray.position, mulConstant(ray.normal, t1));
            if (isIntersecting(intersectionPoint))
            {
                return t1;
            }
        }
        if (t2 > 0)
        {
            Point intersectionPoint = addPoint(ray.position, mulConstant(ray.normal, t2));
            if (isIntersecting(intersectionPoint))
            {
                return t2;
            }
        }

        return -1;
    }
    virtual double helperFunc(Raylight ray, double &color_x, double &color_y, double &color_z, int level)
    {
        auto [C0, C1, C2] = calculateCoefficients(ray);
        double discriminant = calculateDiscriminant(C0, C1, C2);

        if (discriminant < 0 || fabs(C0) < 1e-5)
            return -1;

        auto [t1, t2] = calculateRoots(C0, C1, discriminant);
        if (t1 < 0 && t2 < 0)
            return -1;

        return findIntersection(ray, t1, t2);
    }
};

double determinant(double ara[3][3])
{
    double v1 = ara[0][0] * (ara[1][1] * ara[2][2] - ara[1][2] * ara[2][1]);
    double v2 = ara[0][1] * (ara[1][0] * ara[2][2] - ara[1][2] * ara[2][0]);
    double v3 = ara[0][2] * (ara[1][0] * ara[2][1] - ara[1][1] * ara[2][0]);
    return v1 - v2 + v3;
}

class Triangle : public Material
{
public:
    Point A, B, C;
    Triangle(Point A, Point B, Point C)
    {
        this->A = A;
        this->B = B;
        this->C = C;
    }

    virtual Raylight getNormal(Point point, Raylight incidentRay)
    {
        Point normal = crossProduct(subPoint(B, A), subPoint(C, A));
        normal = normalPoint(normal);

        double dot = mulPoint(normal, incidentRay.normal);

        return Raylight(point, dot < 0 ? mulConstant(normal, -1) : normal);
    }

    virtual void draw()
    {
        glColor3f(color_x, color_y, color_z);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(A.x, A.y, A.z);
            glVertex3f(B.x, B.y, B.z);
            glVertex3f(C.x, C.y, C.z);
        }
        glEnd();
    }
    virtual double helperFunc(Raylight ray, double &color_x, double &color_y, double &color_z, int level)
    {
        double matrices[4][3][3] = {
            {{A.x - ray.position.x, A.x - C.x, ray.normal.x},
             {A.y - ray.position.y, A.y - C.y, ray.normal.y},
             {A.z - ray.position.z, A.z - C.z, ray.normal.z}},

            {{A.x - B.x, A.x - ray.position.x, ray.normal.x},
             {A.y - B.y, A.y - ray.position.y, ray.normal.y},
             {A.z - B.z, A.z - ray.position.z, ray.normal.z}},

            {{A.x - B.x, A.x - C.x, A.x - ray.position.x},
             {A.y - B.y, A.y - C.y, A.y - ray.position.y},
             {A.z - B.z, A.z - C.z, A.z - ray.position.z}},

            {{A.x - B.x, A.x - C.x, ray.normal.x},
             {A.y - B.y, A.y - C.y, ray.normal.y},
             {A.z - B.z, A.z - C.z, ray.normal.z}}};

        double Adet = determinant(matrices[3]);
        double beta = determinant(matrices[0]) / Adet;
        double gamma = determinant(matrices[1]) / Adet;
        double t = determinant(matrices[2]) / Adet;

        return (beta + gamma < 1 && beta > 0 && gamma > 0 && t > 0) ? t : -1;
    }
};

class Sphere : public Material
{
public:
    Sphere(Point center, double radius)
    {
        ref_point = center;
        length = radius;
    }

    virtual Raylight getNormal(Point point, Raylight incidentRay)
    {
        return Raylight(point, subPoint(point, ref_point));
    }

    virtual void draw()
    {
        int stacks = 30;
        int slices = 20;

        Point points[100][100];
        int i, j;
        double h, r;
        for (i = 0; i <= stacks; i++)
        {
            h = length * sin(((double)i / (double)stacks) * (pi / 2));
            r = length * cos(((double)i / (double)stacks) * (pi / 2));
            for (j = 0; j <= slices; j++)
            {
                points[i][j].x = r * cos(((double)j / (double)slices) * 2 * pi);
                points[i][j].y = r * sin(((double)j / (double)slices) * 2 * pi);
                points[i][j].z = h;
            }
        }
        for (i = 0; i < stacks; i++)
        {
            glPushMatrix();
            glTranslatef(ref_point.x, ref_point.y, ref_point.z);
            glColor3f(color_x, color_y, color_z);
            for (j = 0; j < slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);

                    glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
                }
                glEnd();
            }
            glPopMatrix();
        }
    }

    virtual double helperFunc(Raylight ray, double &color_x, double &color_y, double &color_z, int level)
    {

        ray.position = subPoint(ray.position, ref_point);

        double a = 1;
        double b = 2 * mulPoint(ray.normal, ray.position);
        double c = mulPoint(ray.position, ray.position) - (length * length);

        double discriminant = pow(b, 2) - 4 * a * c;
        double t = -1;
        if (discriminant < 0)
        {
            t = -1;
        }
        else
        {

            if (fabs(a) < 1e-5)
            {
                t = -c / b;
                return t;
            }

            double t1 = (-b - sqrt(discriminant)) / (2 * a);
            double t2 = (-b + sqrt(discriminant)) / (2 * a);

            if (t2 < t1)
                swap(t1, t2);

            if (t1 > 0)
            {
                t = t1;
            }
            else if (t2 > 0)
            {
                t = t2;
            }
            else
            {
                t = -1;
            }
        }

        return t;
    }
};

class Floor : public Material
{
public:
    int tiles;

    Floor() : tiles(1) {}

    Floor(int floorWidth, int tileWidth)
    {
        tiles = floorWidth / tileWidth;
        ref_point = Point(-floorWidth / 2, -floorWidth / 2, 0);
        length = tileWidth;
    }

    void getColor(double &color_x, double &color_y, double &color_z, Point point) override
    {
        int tileX = (point.x - ref_point.x) / length;
        int tileY = (point.y - ref_point.y) / length;

        if (tileX < 0 || tileX >= tiles || tileY < 0 || tileY >= tiles || ((tileX + tileY) % 2) != 0)
        {
            color_x = 0;
            color_y = 0;
            color_z = 0;
        }
        else
        {
            color_x = 1;
            color_y = 1;
            color_z = 1;
        }
    }
    Raylight getNormal(Point point, Raylight incidentRay) override
    {
        return Raylight(point, Point(0, 0, incidentRay.normal.z > 0 ? 1 : -1));
    }

    void draw() override
    {
        for (int i = 0; i < tiles; i++)
        {
            for (int j = 0; j < tiles; j++)
            {
                glColor3f(((i + j) % 2) == 0 ? 1 : 0, ((i + j) % 2) == 0 ? 1 : 0, ((i + j) % 2) == 0 ? 1 : 0);

                glBegin(GL_QUADS);
                {
                    glVertex3f(ref_point.x + i * length, ref_point.y + j * length, 0);
                    glVertex3f(ref_point.x + (i + 1) * length, ref_point.y + j * length, 0);
                    glVertex3f(ref_point.x + (i + 1) * length, ref_point.y + (j + 1) * length, 0);
                    glVertex3f(ref_point.x + i * length, ref_point.y + (j + 1) * length, 0);
                }
                glEnd();
            }
        }
    }

     double helperFunc(Raylight ray, double &color_x, double &color_y, double &color_z, int level) override
    {
        Point normal = Point(0, 0, 1);
        double dotP = mulPoint(normal, ray.normal);

        if (round(dotP * 100) == 0)
            return -1;

        double t = -mulPoint(normal, ray.position) / dotP;

        Point p = addPoint(ray.position, mulConstant(ray.normal, t));

        if (p.x <= ref_point.x || p.x >= abs(ref_point.x) || p.y <= ref_point.y || p.y >= abs(ref_point.y))
        {
            return -1;
        }

        return t;
    }
};