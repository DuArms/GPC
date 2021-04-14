#include <GL/glut.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream>

#define or ||
#define and &&

#define x first
#define y second

// dimensiunea ferestrei in pixeli
#define dim 1000

unsigned char prevKey;

enum EObiect {cubw, cubs, sferaw, sferas, triangle, wireCube} ob = cubw;

double EPS = 1e-6;
double PI = acos(-1);

double radiansToAngle(double radians){
    return radians / (2*PI) * 360;
}

double angleToRadians(double angle){
    return angle/360*2*PI;
}

void crossProduct(double vect_A[], double vect_B[], double cross_P[])
{

    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}

double dotProduct(double vect_A[], double vect_B[])
{
    double product = 0.0;

    for (int i = 0; i < 3; i++)

        product = product + vect_A[i] * vect_B[i];
    return product;
}

void P2DToPolar(double x, double y, double &r, double &t) {
    if (abs(x) < EPS and abs(y) < EPS){
        r = t = 0.0;
        return;
    }

    r = sqrt(x*x + y*y);
    t = atan2(y, x);
}

void P2DRotAxis(double &x, double &y, double alpha){
    double t, r;
    P2DToPolar(x, y, r, t);

    x = r*cos(t + alpha);
    y = r*sin(t + alpha);
}

struct Point{
    double x{};
    double y{};
    double z{};

    double r{};
    double theta{}; // rotation along Z
    double rho{}; //

    Point(double ax, double ay, double az) : x(ax), y(ay), z(az) {CartesianToPolar();}
    Point() : x(0.0), y(0.0), z(0.0), r(0.0), theta(0.0), rho(0.0) {}

    void CartesianToPolar() {
        if (abs(x) < EPS and abs(y) < EPS and abs(z) < EPS) {
            r = theta = rho = 0.0;
            return;
        }

        r = sqrt(x*x + y*y + z*z);
        theta = acos(z / r);
        rho = atan2(y, x);
    }

    void PolarToCartesian() {
        x = r * cos(rho) * sin(theta);
        y = r * sin(rho) * sin(theta);
        z = r*cos(theta);
    }

    void debug(const std::string& name="P") const {
        std::cout << name << ": [";
        std::cout << "x=" << x << ", ";
        std::cout << "y=" << y << ", ";
        std::cout << "z=" << z << " | ";

        std::cout << "r=" << r << ", ";
        std::cout << "a=" << radiansToAngle(theta) << ", ";
        std::cout << "b=" << radiansToAngle(rho) << "]\n";
    }

    Point operator- (Point const& other){
        Point res;
        res.x = x - other.x;
        res.y = y - other.y;
        res.z = z - other.z;
        return res;
    }
};

void DisplayAxe(double offSet = 0.0) {
    int X, Y, Z;
    X = Y = 200;
    Z = 200;

    glLineWidth(2);

    // axa Ox - verde
    glColor3f(0, 1, 0);
    glBegin(GL_LINE_STRIP);
    glVertex3f(offSet,offSet,offSet);
    glVertex3f(X,0,0);
    glEnd();

    // axa Oy - albastru
    glColor3f(0, 0, 1);
    glBegin(GL_LINE_STRIP);
    glVertex3f(offSet,offSet,offSet);
    glVertex3f(0,Y,0);
    glEnd();

    // axa Oz - rosu
    glColor3f(1, 0, 0);
    glBegin(GL_LINE_STRIP);
    glVertex3f(offSet,offSet,offSet);
    glVertex3f(0,0,Z);
    glEnd();

    glColor3f(1, 0, 1);
    glBegin(GL_LINE_STRIP);
    glVertex3f(offSet,offSet,offSet);
    glVertex3f(Z,Z,Z);
    glEnd();

    glLineWidth(1);
}

// cub wireframe
void Display1() {
    glColor3f(1,0,0);
    glutWireCube(1);
}

// cub solid
void Display2() {
    glColor3f(1,0,0);
    glutSolidCube(1);
}

// sfera wireframe
void Display3() {
    glColor3f(0,0,1);
    glutWireSphere(1, 10, 10);
}

// sfera solida
void Display4() {
    glColor3f(0,0,1);
    glutSolidSphere(1, 10, 10);
}

void Display5(bool rot = false) {
    glColor3f(0,0,1);

    Point p[3];

    p[0] = Point(1.0, 0.0, 1.0);
    p[1] = Point(2.0, 0.0, 2.0);
    p[2] = Point(1.0, 1.0, 1.0);

    p[0] = Point(1.0, 2.0, 3.0);
    p[1] = Point(3.0, 1.0, 2.0);
    p[2] = Point(8.0, 8.0, -8.0);

    p[0].debug("P0");
    p[1].debug("P1");
    p[2].debug("P2");

    glPushMatrix();
    {
        Point pRot[3];
        for (int i = 0; i < 3; i++) {
            pRot[i] = p[i] - p[0];

            std::stringstream ss;
            ss << i;
            std::string s;
            ss >> s;
            pRot[i].debug("PRot Trans" + s);
        }
        std::cout << "\n";

        double r, t;

        double r1, t1;
        P2DToPolar(pRot[1].z, pRot[1].x, r1, t1);

        for (int i = 0; i < 3; i++){
            P2DRotAxis(pRot[i].z, pRot[i].x, -t1);

            std::stringstream ss;
            ss << i;
            std::string s;
            ss >> s;
            pRot[i].debug("PRot ZX" + s);
        }
        std::cout << "\n";

        double r2, t2;
        P2DToPolar(pRot[1].y, pRot[1].z, r2, t2);

        for (int i = 0; i < 3; i++){
            P2DRotAxis(pRot[i].y, pRot[i].z, -t2);
            P2DRotAxis(pRot[i].y, pRot[i].z, PI/2);

            std::stringstream ss;
            ss << i;
            std::string s;
            ss >> s;
            pRot[i].debug("PRot YZ" + s);
        }
        std::cout << "\n";

        double r3, t3;
        P2DToPolar(pRot[2].x, pRot[2].y, r3, t3);
        for (int i = 0; i < 3; i++){
            P2DRotAxis(pRot[i].x, pRot[i].y, -t3);
            P2DRotAxis(pRot[i].x, pRot[i].y, PI/2);

            std::stringstream ss;
            ss << i;
            std::string s;
            ss >> s;
            pRot[i].debug("PRot XY" + s);
        }
        std::cout << "\n";

        glRotated(radiansToAngle(PI/2), 0.0, 0.0, r3);
        glRotated(-radiansToAngle(t3), 0.0, 0.0, r3);

        glRotated(radiansToAngle(PI/2), r2, 0.0, 0.0);
        glRotated(-radiansToAngle(t2), r2, 0.0, 0.0);

        glRotated(-radiansToAngle(t1), 0.0, r1, 0.0);
        glTranslatef(-p[0].x, -p[0].y, -p[0].z);

        glBegin(GL_TRIANGLE_FAN);
        for (auto &i : p) {
            glVertex3f(i.x, i.y, i.z);
        }
        glEnd();
    }
    glPopMatrix();
}

double angle = 0.0;

void Display6(){
    glColor3f(0,0,1);

    Point pts[8];
    for (uint64_t i = 0; i < 8; i++){
        double x = i & 1U;
        double y = (i & 2U) >> 1;
        double z = (i & 4U) >> 2;
        pts[i] = Point(x, y, z);
    }
    std::cout << "\n";

    std::vector<Point> faces[6];
    faces[0] = {pts[0], pts[1], pts[3], pts[2]};
    faces[1] = {pts[0], pts[1], pts[5], pts[4]};
    faces[2] = {pts[1], pts[3], pts[7], pts[5]};
    faces[3] = {pts[0], pts[2], pts[6], pts[4]};
    faces[4] = {pts[2], pts[3], pts[7], pts[6]};
    faces[5] = {pts[4], pts[5], pts[7], pts[6]};

    glPushMatrix();
    {
        double r1, t1;
        Point rotPt = pts[7];
        P2DToPolar(rotPt.z, rotPt.x, r1, t1);
        P2DRotAxis(rotPt.z, rotPt.x, -t1);
        rotPt.debug("PRotZX");

        double r2, t2;
        P2DToPolar(rotPt.y, rotPt.z, r2, t2);
        P2DRotAxis(rotPt.y, rotPt.z, -t2);
        rotPt.debug("PRotYZ");

        //Rotations
        glRotated(radiansToAngle(t1), 0.0, r1, 0.0);
        glRotated(radiansToAngle(t2), r2, 0.0, 0.0);

        glRotated(angle, 0.0, 1.0, 0.0);

        glRotated(-radiansToAngle(t2), r2, 0.0, 0.0);
        glRotated(-radiansToAngle(t1), 0.0, r1, 0.0);

        glBegin(GL_LINES);
        for (auto const& face : faces) {
            for (int i = 0; i < 4; i++) {
                int j = i < 3 ? i + 1 : 0;
                glVertex3f(face[i].x, face[i].y, face[i].z);
                glVertex3f(face[j].x, face[j].y, face[j].z);
            }
        }
        glEnd();
    }
    glPopMatrix();
}

void DisplayObiect()
{
    switch (ob)
    {
        case cubw:
            Display1();
            break;

        case cubs:
            Display2();
            break;

        case sferaw:
            Display3();
            break;

        case sferas:
            Display4();
            break;

        case triangle:
            Display5();
            break;

        case wireCube:
            Display6();
            break;

        default:
            break;
    }
}

// rotatia cu un unghi de 10 grade in raport cu axa x
void DisplayX() {
    glMatrixMode(GL_MODELVIEW);
    glRotated(10,1,0,0);
}

// rotatia cu un unghi de 10 grade in raport cu axa y
void DisplayY() {
    glMatrixMode(GL_MODELVIEW);
    glRotated(10,0,1,0);
}

// rotatia cu un unghi de 10 grade in raport cu axa z
void DisplayZ() {
    glMatrixMode(GL_MODELVIEW);
    glRotated(10,0,0,1);
}

// Translatia cu 0.2, 0.2, 0.2
void DisplayT() {
    glMatrixMode(GL_MODELVIEW);
    glTranslatef(0.2, 0.2, 0.2);
}

// Scalarea cu 1.2, 1.2, 1.2
void DisplayS() {
    glMatrixMode(GL_MODELVIEW);
    glScalef(1.2, 1.2, 1.2);
}

void Init() {
    glClearColor(1, 1, 1, 1);
    glLineWidth(2);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-10, 10, -10, 10, 30, -30);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glRotated(20, 1, 0, 0);
    glRotated(-20, 0, 1, 0);
}

void Display() {
    switch(prevKey)
    {
        case 'a':
            DisplayAxe();
            break;

        case '0':
            glClear(GL_COLOR_BUFFER_BIT);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            glRotated(20, 1, 0, 0);
            glRotated(-20, 0, 1, 0);
            break;

        case '1':
            Display1();
            ob = cubw;
            break;

        case '2':
            Display2();
            ob = cubs;
            break;

        case '3':
            Display3();
            ob = sferaw;
            break;

        case '4':
            Display4();
            ob = sferas;
            break;

        case '5':
            Display5();
            ob = triangle;
            break;

        case '6': {
            glClear(GL_COLOR_BUFFER_BIT);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            double l = 5.0;
            gluLookAt(l, l, l,
                      0.0, 0.0, 0.0,
                      0.0, l, 0.0
            );
//            gluLookAt(l, 0.0, l,
//                      0.0, 0.0, 0.0,
//                      0.0, l, 0.0
//            );
//            glRotated(45, 1, 0, 0);
//            glRotated(45, 0, 0, 1);
            ob = wireCube;
            Display6();
        }
            break;

        case 'q':
            glClear(GL_COLOR_BUFFER_BIT);
            angle += 10.0;
            DisplayAxe();
            DisplayObiect();
            break;

        case 'x':
            glClear(GL_COLOR_BUFFER_BIT);
            DisplayX();
            DisplayAxe();
            DisplayObiect();
            break;

        case 'y':
            glClear(GL_COLOR_BUFFER_BIT);
            DisplayY();
            DisplayAxe();
            DisplayObiect();
            break;

        case 'z':
            glClear(GL_COLOR_BUFFER_BIT);
            DisplayZ();
            DisplayAxe();
            DisplayObiect();
            break;

        case 't':
            glClear(GL_COLOR_BUFFER_BIT);
            DisplayT();
            DisplayAxe();
            DisplayObiect();
            break;

        case 's':
            glClear(GL_COLOR_BUFFER_BIT);
            DisplayS();
            DisplayAxe();
            DisplayObiect();
            break;

        default:
            break;
    }
    glutSwapBuffers();
}

void Reshape(int w, int h) {
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}

void KeyboardFunc(unsigned char key, int x, int y) {
    prevKey = key;
    if (key == 27) // escape
        exit(0);
    glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y) {
}

int main(int argc, char** argv) {

    glutInit(&argc, argv);

    glutInitWindowSize(dim, dim);

    glutInitWindowPosition(100, 100);

    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);

    glutCreateWindow (argv[0]);

    Init();

    glutReshapeFunc(Reshape);

    glutKeyboardFunc(KeyboardFunc);

    glutMouseFunc(MouseFunc);

    glutDisplayFunc(Display);

    glutMainLoop();

    return 0;
}
