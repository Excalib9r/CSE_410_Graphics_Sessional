#include <bits/stdc++.h>
#include <math.h>
#include "bitmap_image.hpp"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Point
{
public:
    double x, y, z;
    Point()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    Point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
};

static unsigned long int g_seed = 1;
inline int rand()
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

void updateIntersectPoint1(Point &p, double y, pair<Point, Point> &intersect_points)
{
    double minX = min(p.x, intersect_points.first.x);
    double maxX = max(p.x, intersect_points.second.x);

    intersect_points.first = (minX == p.x) ? Point(p.x, y, p.z) : intersect_points.first;
    intersect_points.second = (maxX == p.x) ? Point(p.x, y, p.z) : intersect_points.second;
}

void updateIntersectPoint2(Point& p1, Point& p2, double y, pair<Point, Point>& intersect_points) {
    double x = p1.x;
    double z = p1.z;

    if (p1.y != p2.y)
    {
        x = p1.x + (y - p1.y) * (p1.x - p2.x) / (p1.y - p2.y);
        z = p1.z + (y - p1.y) * (p1.z - p2.z) / (p1.y - p2.y);
    }

    double minX = std::min(p1.x, p2.x);
    double maxX = std::max(p1.x, p2.x);

    if (x >= minX && x <= maxX)
    {
        intersect_points.first = (x < intersect_points.first.x) ? Point(x, y, z) : intersect_points.first;
        intersect_points.second = (x > intersect_points.second.x) ? Point(x, y, z) : intersect_points.second;
    }
}

void getIntersectingPoints(Point p1, Point p2, pair<Point, Point> &intersect_points, double y)
{
    if (abs(p1.y - p2.y) < 1e-9)
    {
        if (abs(p1.y - y) < 1e-9)
        {
            updateIntersectPoint1(p1, y, intersect_points);
            updateIntersectPoint1(p2, y, intersect_points);
        }
    }
    else
    {
        updateIntersectPoint2(p1, p2, y, intersect_points);
    }
}

class Triangle
{
public:
    Point p1, p2, p3;
    int color[3];
    Triangle()
    {
        this->p1 = Point();
        this->p2 = Point();
        this->p3 = Point();

        for (int i = 0; i < 3; i++)
        {
            color[i] = rand() % 256;
        }
    }
    Triangle(Point p1, Point p2, Point p3)
    {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;

        for (int i = 0; i < 3; i++)
        {
            color[i] = rand() % 256;
        }
    }
    void setPoints(Point p1, Point p2, Point p3)
    {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
    }
    void printTriangle(std::ofstream &file)
    {
        file << std::fixed << std::setprecision(8) << p1.x << " " << p1.y << " " << p1.z << endl;
        file << std::fixed << std::setprecision(8) << p2.x << " " << p2.y << " " << p2.z << endl;
        file << std::fixed << std::setprecision(8) << p3.x << " " << p3.y << " " << p3.z << endl;
        file << endl;
    }

    pair<Point, Point> findIntersectingPoints(double y)
    {
        pair<Point, Point> intersect_points;
        intersect_points.first = Point(1e18, y, 1e18);
        intersect_points.second = Point(-1e18, y, -1e18);
        getIntersectingPoints(p1, p2, intersect_points, y);
        getIntersectingPoints(p2, p3, intersect_points, y);
        getIntersectingPoints(p3, p1, intersect_points, y);
        return intersect_points;
    }
};

void Normalize(double *x, double *y, double *z)
{
    double length = sqrt((*x) * (*x) + (*y) * (*y) + (*z) * (*z));
    *x = (*x) / length;
    *y = (*y) / length;
    *z = (*z) / length;
}

class Rotate
{
public:
    double angle, x, y, z;
    Rotate(double angle, double x, double y, double z)
    {
        this->angle = angle;
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Rotate()
    {
        this->angle = 0;
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    void setRotate(double angle, double x, double y, double z)
    {
        this->angle = angle;
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void printRotate()
    {
        std::cout << angle << " " << x << " " << y << " " << z << endl;
    }
};

class Translate
{
public:
    double x, y, z;
    Translate(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Translate()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    void setTranslate(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void printTranslate()
    {
        std::cout << x << " " << y << " " << z << endl;
    }
};

class Scale
{
public:
    double x, y, z;
    Scale(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Scale()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    void setScale(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void printScale()
    {
        std::cout << x << " " << y << " " << z << endl;
    }
};

void traingleToMatraix(Triangle *t, double **mat)
{
    mat[0][0] = t->p1.x;
    mat[0][1] = t->p2.x;
    mat[0][2] = t->p3.x;
    mat[1][0] = t->p1.y;
    mat[1][1] = t->p2.y;
    mat[1][2] = t->p3.y;
    mat[2][0] = t->p1.z;
    mat[2][1] = t->p2.z;
    mat[2][2] = t->p3.z;
    mat[3][0] = 1;
    mat[3][1] = 1;
    mat[3][2] = 1;
}

void matrixToTriangle(double **mat, Triangle *t)
{
    t->p1.x = mat[0][0];
    t->p2.x = mat[0][1];
    t->p3.x = mat[0][2];
    t->p1.y = mat[1][0];
    t->p2.y = mat[1][1];
    t->p3.y = mat[1][2];
    t->p1.z = mat[2][0];
    t->p2.z = mat[2][1];
    t->p3.z = mat[2][2];
}

class Matrix
{
public:
    double mat[4][4];
    Matrix(double mat[4][4])
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                this->mat[i][j] = mat[i][j];
            }
        }
    }
    Matrix(Translate t)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                if (i == j)
                    mat[i][j] = 1;
                else if (i == 0 && j == 3)
                    mat[i][j] = t.x;
                else if (i == 1 && j == 3)
                    mat[i][j] = t.y;
                else if (i == 2 && j == 3)
                    mat[i][j] = t.z;
                else
                    mat[i][j] = 0;
            }
        }
    }
    Matrix(Rotate r)
    {
        double angle = r.angle * M_PI / 180.0;
        ;
        double x = r.x;
        double y = r.y;
        double z = r.z;
        double c = cos(angle);
        if (std::abs(c) < 1e-15)
        {
            c = 0;
        }
        double s = sin(angle);
        if (std::abs(s) < 1e-15)
        {
            s = 0;
        }
        double t = 1 - c;
        Normalize(&x, &y, &z);
        double arr[3][3] = {{x * x * t + c, x * y * t - z * s, x * z * t + y * s},
                            {y * x * t + z * s, y * y * t + c, y * z * t - x * s},
                            {z * x * t - y * s, z * y * t + x * s, z * z * t + c}};
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                mat[i][j] = arr[i][j];
            }
        }
        for (int i = 0; i < 3; i++)
        {
            mat[i][3] = 0;
            mat[3][i] = 0;
        }
        mat[3][3] = 1;
    }
    Matrix(Scale s)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                if (i == 0 && i == j)
                    mat[i][j] = s.x;
                else if (i == 1 && i == j)
                    mat[i][j] = s.y;
                else if (i == 2 && i == j)
                    mat[i][j] = s.z;
                else if (i == j)
                    mat[i][j] = 1;
                else
                    mat[i][j] = 0;
            }
        }
    }
    void matrixMul(Matrix *m1)
    {
        double dummy[4][4];
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                dummy[i][j] = 0;
                for (int k = 0; k < 4; k++)
                {
                    dummy[i][j] += m1->mat[i][k] * mat[k][j];
                }
            }
        }
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                m1->mat[i][j] = dummy[i][j];
            }
        }
    }
    void printMatrix()
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                std::cout << mat[i][j] << " ";
            }
            std::cout << endl;
        }
    }
};

void MatrixMultiplication(Matrix *mat1, double **mat2, double **mat3)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            mat3[i][j] = 0;
            for (int k = 0; k < 4; k++)
            {
                mat3[i][j] += mat1->mat[i][k] * mat2[k][j];
            }
        }
    }
}

void normalizeResult(double **mat)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            mat[i][j] = mat[i][j] / mat[3][j];
        }
    }
}

double get_z_value(double a, double b, double c, double d, double x)
{
    return (b != a)
               ? c + (x - a) * (d - c) / (b - a)
               : c;
}

void processPixel(double x, double y, double del_x, double del_y, pair<Point, Point>& interSect, vector<vector<double>>& z_buffer, vector<pair<int, int>>& f_buffer) {
    double z = get_z_value(interSect.first.x, interSect.second.x, interSect.first.z, interSect.second.z, x);
    int i = ((x + 1.0) / del_x);
    int j = ((1.0 - y) / del_y);
    if (z + 1.0 > 1e-9 && z_buffer[i][j] - z > 1e-9)
    {
        z_buffer[i][j] = z;
        f_buffer.push_back(make_pair(i, j));
    }
}

double get_topScanline(Triangle *t, double upper, double del_y)
{
    double topScanline = max({t->p1.y, t->p2.y, t->p3.y});
    int rows = round((1 - topScanline) / del_y);
    topScanline = (upper - topScanline < 1e-9) ? upper : upper - rows * del_y;
    return topScanline;
}

double get_bottomScanline(Triangle *t, double lower, double del_y)
{
    double bottomScanline = min({t->p1.y, t->p2.y, t->p3.y});
    int rows = round((bottomScanline + 1) / del_y);
    bottomScanline = (bottomScanline - lower < 1e-9) ? lower : lower + rows * del_y;
    return bottomScanline;
}


void processIntersectingPoints(Triangle* t, double y, double del_x, double del_y, double leftX, double rightX, vector<vector<double>>& z_buffer, vector<pair<int, int>>& f_buffer) {
    pair<Point, Point> interSect = t->findIntersectingPoints(y);

    double left_interSect_col = interSect.first.x;
    double right_interSect_col = interSect.second.x;

    int cols = ((left_interSect_col + 1.0) / del_x);
    left_interSect_col = (left_interSect_col - leftX < 1e-9) ? leftX : leftX + cols * del_x;

    cols = ((1.0 - right_interSect_col) / del_x);
    right_interSect_col = (rightX - right_interSect_col < 1e-9) ? rightX : rightX - cols * del_x;

    for (double x = left_interSect_col; (x - right_interSect_col) < 1e-9; x += del_x)
    {   
        processPixel(x, y, del_x, del_y, interSect, z_buffer, f_buffer);
    }
}

double eye_x, eye_y, eye_z;
double look_x, look_y, look_z;
double up_x, up_y, up_z;
double fovY, aspect_ratio, near, far;

int main()
{
    std::ifstream file("5/scene.txt");
    std::ifstream config("5/config.txt");
    std::ofstream file2("out/stage1.txt");
    std::ofstream file3("out/stage2.txt");
    std::ofstream file4("out/stage3.txt");
    std::ofstream file5("out/z_buffer.txt");
    std::string line;
    int counter = 0;
    std::stack<Matrix *> matrices;
    std::vector<Triangle *> triangles;
    double identity[4][4] = {{1, 0, 0, 0},
                             {0, 1, 0, 0},
                             {0, 0, 1, 0},
                             {0, 0, 0, 1}};
    matrices.push(new Matrix(identity));
    Matrix *transform = matrices.top();
    double l_x, l_y, l_z;
    double r_x, r_y, r_z;
    double u_x, u_y, u_z;

    Matrix *projection = new Matrix(identity);
    Matrix *T = new Matrix(identity);
    Matrix *view = new Matrix(identity);

    while (std::getline(file, line))
    {
        if (counter < 3)
        {
            std::istringstream iss(line);
            string a, b, c;
            if (!(iss >> a >> b >> c))
            {
                break;
            } // error
            if (counter == 0)
            {
                eye_x = stod(a);
                eye_y = stod(b);
                eye_z = stod(c);
                double t[4][4] = {{1, 0, 0, -eye_x},
                                  {0, 1, 0, -eye_y},
                                  {0, 0, 1, -eye_z},
                                  {0, 0, 0, 1}};
                delete T;
                T = new Matrix(t);
            }
            else if (counter == 1)
            {
                look_x = stod(a);
                look_y = stod(b);
                look_z = stod(c);
            }
            else
            {
                up_x = stod(a);
                up_y = stod(b);
                up_z = stod(c);

                l_x = look_x - eye_x;
                l_y = look_y - eye_y;
                l_z = look_z - eye_z;
                Normalize(&l_x, &l_y, &l_z);
                r_x = l_y * up_z - l_z * up_y;
                r_y = l_z * up_x - l_x * up_z;
                r_z = l_x * up_y - l_y * up_x;
                Normalize(&r_x, &r_y, &r_z);
                u_x = r_y * l_z - r_z * l_y;
                u_y = r_z * l_x - r_x * l_z;
                u_z = r_x * l_y - r_y * l_x;
                Normalize(&u_x, &u_y, &u_z);

                double V[4][4] = {{r_x, r_y, r_z, 0},
                                  {u_x, u_y, u_z, 0},
                                  {-l_x, -l_y, -l_z, 0},
                                  {0, 0, 0, 1}};
                delete view;
                view = new Matrix(V);
                T->matrixMul(view);
            }
        }
        else if (counter == 3)
        {
            std::istringstream iss(line);
            string a, b, c, d;
            if (!(iss >> a >> b >> c >> d))
            {
                break;
            } // error
            fovY = stod(a);
            aspect_ratio = stod(b);
            near = stod(c);
            far = stod(d);
            double fovX = fovY * aspect_ratio;
            double t = near * tan(fovY / 2 * M_PI / 180.0);
            double r = near * tan(fovX / 2 * M_PI / 180.0);
            double p[4][4] = {{near / r, 0, 0, 0},
                              {0, near / t, 0, 0},
                              {0, 0, -(far + near) / (far - near), -(2 * far * near) / (far - near)},
                              {0, 0, -1, 0}};
            delete projection;
            projection = new Matrix(p);
        }
        else
        {
            std::istringstream iss(line);
            string val;
            if (!(iss >> val))
            {
                break;
            } // error
            if (val == "triangle")
            {
                std::getline(file, line);
                std::istringstream iss2(line);
                string a1, a2, a3;
                if (!(iss2 >> a1 >> a2 >> a3))
                {
                    break;
                } // error
                std::getline(file, line);
                std::istringstream iss3(line);
                string b1, b2, b3;
                if (!(iss3 >> b1 >> b2 >> b3))
                {
                    break;
                } // error
                std::getline(file, line);
                std::istringstream iss4(line);
                string c1, c2, c3;
                if (!(iss4 >> c1 >> c2 >> c3))
                {
                    break;
                } // error
                Triangle *t = new Triangle();
                Point p1 = Point(stod(a1), stod(a2), stod(a3));
                Point p2 = Point(stod(b1), stod(b2), stod(b3));
                Point p3 = Point(stod(c1), stod(c2), stod(c3));
                t->setPoints(p1, p2, p3);
                double **mat2 = new double *[4];
                for (int i = 0; i < 4; i++)
                {
                    mat2[i] = new double[3];
                }
                double **result = new double *[4];
                for (int i = 0; i < 4; i++)
                {
                    result[i] = new double[3];
                }

                traingleToMatraix(t, mat2);
                MatrixMultiplication(transform, mat2, result);
                normalizeResult(result);
                matrixToTriangle(result, t);
                t->printTriangle(file2);
                traingleToMatraix(t, mat2);
                MatrixMultiplication(view, mat2, result);
                normalizeResult(result);
                matrixToTriangle(result, t);
                t->printTriangle(file3);
                traingleToMatraix(t, mat2);
                MatrixMultiplication(projection, mat2, result);
                normalizeResult(result);
                matrixToTriangle(result, t);
                t->printTriangle(file4);
                triangles.push_back(t);
    
                for(int i = 0; i < 4; i++){
                    delete[] mat2[i];
                    delete[] result[i];
                }
                delete[] mat2;
                delete[] result;
            }
            else if (val == "translate")
            {
                std::getline(file, line);
                std::istringstream iss2(line);
                string a1, a2, a3;
                if (!(iss2 >> a1 >> a2 >> a3))
                {
                    break;
                } // error
                Translate t = Translate();
                t.setTranslate(stod(a1), stod(a2), stod(a3));
                Matrix m = Matrix(t);
                m.matrixMul(transform);
            }
            else if (val == "rotate")
            {
                std::getline(file, line);
                std::istringstream iss2(line);
                string a1, a2, a3, a4;
                if (!(iss2 >> a1 >> a2 >> a3 >> a4))
                {
                    break;
                } // error
                Rotate r = Rotate();
                r.setRotate(stod(a1), stod(a2), stod(a3), stod(a4));
                Matrix m = Matrix(r);
                m.matrixMul(transform);
            }
            else if (val == "scale")
            {
                std::getline(file, line);
                std::istringstream iss2(line);
                string a1, a2, a3;
                if (!(iss2 >> a1 >> a2 >> a3))
                {
                    break;
                } // error
                Scale s = Scale();
                s.setScale(stod(a1), stod(a2), stod(a3));
                Matrix m = Matrix(s);
                m.matrixMul(transform);
            }
            else if (val == "push")
            {
                Matrix *m = new Matrix(*transform);
                matrices.push(m);
            }
            else if (val == "pop")
            {
                transform = matrices.top();
                matrices.pop();
            }
            else
            {
                break;
            }
        }
        counter++;
    }

    delete projection;
    delete T;
    delete view;
    delete transform;

    file.close();
    file2.close();
    file3.close();
    file4.close();

    string screen_width, screen_height;
    double del_x, del_y, upper, lower, leftX, rightX;

    config >> screen_width >> screen_height;
    double s_height = stod(screen_height);
    double s_width = stod(screen_width);
    bitmap_image image = bitmap_image(s_width, s_height);

    vector<vector<double>> z_buffer(s_height, vector<double>(s_width, 1));
    vector<vector<pair<int, int>>> f_buffer(triangles.size());

    del_x = 2.0 / s_width;
    del_y = 2.0 / s_height;

    upper = 1.0 - del_y / 2.0;
    lower = -1.0 + del_y / 2.0;
    leftX = -1.0 + del_x / 2.0;
    rightX = 1.0 - del_x / 2.0;

    int counter2 = 0;
    for (Triangle *t : triangles)
    {   
        vector<pair<int, int>> helper;
        f_buffer.push_back(helper);

        double top_scanline = get_topScanline(t, upper, del_y);
        double bottom_scanline = get_bottomScanline(t, lower, del_y);

        for (double y = bottom_scanline; (y - top_scanline) < 1e-9; y += del_y)
        {
            processIntersectingPoints(t, y, del_x, del_y, leftX, rightX, z_buffer, f_buffer[counter2]);
        }
        counter2++;
    }

    for (int i = 0; i < f_buffer.size(); i++)
    {
        for (int j = 0; j < f_buffer[i].size(); j++)
        {
            int x = f_buffer[i][j].first;
            int y = f_buffer[i][j].second;
            image.set_pixel(x, y, triangles[i]->color[0], triangles[i]->color[1], triangles[i]->color[2]);
        }
    }

    for (int j = 0; j < s_height; j++)
    {
        for (int i = 0; i < s_width; i++)
        {
            if (abs(z_buffer[i][j] - 1) > 1e-9)
                file5 << z_buffer[i][j] << "\t";
        }
        file5 << endl;
    }

    image.save_image("out/out.bmp");
    for(int i = 0; i < triangles.size(); i++){
        delete triangles[i];
    }

    file5.close();
    config.close();

    std::cout << "Done" << endl;
    return 0;
}
