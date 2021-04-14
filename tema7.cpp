#include <GL/glut.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
// dimensiunea ferestrei in pixeli
#define dim 300

#define or ||
#define  and &&
double EPS = 1e-6;
double PI = acos(-1);

unsigned char prevKey;
GLint k;
// latura cubului
GLdouble lat = 5;

void proiectieParalela(unsigned char);

void DisplayAxe();

void InitObiect();

void DisplayObiect();

double radiansToAngle(double radians) {
    return radians / (2 * PI) * 360;
}

double angleToRadians(double angle) {
    return angle / 360 * 2 * PI;
}

void crossProduct(double vect_A[], double vect_B[], double cross_P[]) {

    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}

double dotProduct(double vect_A[], double vect_B[]) {
    double product = 0.0;

    for (int i = 0; i < 3; i++)

        product = product + vect_A[i] * vect_B[i];
    return product;
}

void P2DToPolar(double x, double y, double &r, double &t) {
    if (abs(x) < EPS and abs(y) < EPS) {
        r = t = 0.0;
        return;
    }

    r = sqrt(x * x + y * y);
    t = atan2(y, x);
}

void P2DRotAxis(double &x, double &y, double alpha) {
    double t, r;
    P2DToPolar(x, y, r, t);

    x = r * cos(t + alpha);
    y = r * sin(t + alpha);
}

void Init(void) {
    glClearColor(1, 1, 1, 1);

    // validare test de adancime
    glEnable(GL_DEPTH_TEST);

    // se aloca 1 lista de display numerotata k
    k = glGenLists(1);
    InitObiect();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
}

void Img1() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1, 7, -7, +1, -5, 2);
    glRotatef(90, 1, 0, 0);

    glMatrixMode(GL_MODELVIEW);
    DisplayAxe();
    DisplayObiect();
}

void Img2() {
    Init();
    glPushMatrix();
    {
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        {
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();

            double l = 1;
            glFrustum(-l, l, -l, l, 1, 100);

            glMatrixMode(GL_MODELVIEW);
            glTranslated(-7.5 - 4, -7.5 - 4, -20);

            glPushMatrix();
            {
                DisplayAxe();
                DisplayObiect();
            }
            glPopMatrix();
        }
        glPopMatrix();
    }
    glPopMatrix();

}


void Img3() {
    Init();
    glPushMatrix();
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        double alfa = 135 + 90 - 15;
        GLdouble cavalier[] =
                {
                        1, 0, 0, 0,
                        0, 1, 0, 0,
                        -1 * cos(angleToRadians(alfa)), -1 * sin(angleToRadians(alfa)), 1, 0,
                        0, 0, 0, 1
                };
        double ratio = 1;
        double w = 30;

        // glLoadMatrixd(cavalier);
        glMultMatrixd(cavalier);
        glOrtho(-w * ratio, w * ratio, -w, w, -30, 100);


        glScaled(3,3,3);
        gluLookAt(2.5, 2.5, 5,
                  2.5, 2.5, 0,
                  0, 1, 0);

        DisplayAxe();
        DisplayObiect();


    }
    glPopMatrix();
};


void Display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    switch (prevKey) {
        case '1':
            Img1();
            break;
        case '2':
            Img2();
            break;
        case '3':
            Img3();
            break;
        case '0':
            // resetarea stivei matricilor de modelare si a celei de proiectie
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            DisplayAxe();
            break;
        case 'a':
            DisplayAxe();
            break;
        case 'c':
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            break;
        case 'x':
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            proiectieParalela('X');
            glMatrixMode(GL_MODELVIEW);
            glRotatef(10, 1, 0, 0);
            DisplayAxe();
            DisplayObiect();
            break;
        case 'y':
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            proiectieParalela('Y');
            glMatrixMode(GL_MODELVIEW);
            glRotatef(10, 0, 1, 0);
            DisplayAxe();
            DisplayObiect();
            break;
        case 'z':
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            proiectieParalela('Z');
            glMatrixMode(GL_MODELVIEW);
            glRotatef(10, 0, 0, 1);
            DisplayAxe();
            DisplayObiect();
            break;
        case 'q':
            // proiectie paralela ortografica frontala (fata)
            proiectieParalela('q');
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            DisplayAxe();
            glTranslated(0, 0, -lat);
            DisplayObiect();
            break;
        case 'w':
            // proiectie paralela ortografica frontala (spate)
            proiectieParalela('w');
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            glPushMatrix();
            glTranslated(0, 0, -lat);
            glTranslated(lat / 2.0, lat / 2.0, lat / 2.0);
            glRotated(180, 0, 1, 0);
            glTranslated(-lat / 2.0, -lat / 2.0, -lat / 2.0);
            DisplayAxe();
            DisplayObiect();
            glPopMatrix();
            break;
        default:
            break;
    }
    glutSwapBuffers();

}

void Reshape(int w, int h) {
    h = (h == 0) ? 1 : h;
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

int main(int argc, char **argv) {

    glutInit(&argc, argv);

    glutInitWindowSize(dim, dim);

    glutInitWindowPosition(100, 100);

    // GLUT_DEPTH - eliminarea suprafetelor ascunse cu alg. Z-buffer
    glutInitDisplayMode(GL_COLOR_BUFFER_BIT | GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    glutCreateWindow(argv[0]);

    Init();

    glutReshapeFunc(Reshape);

    glutKeyboardFunc(KeyboardFunc);

    glutMouseFunc(MouseFunc);

    glutDisplayFunc(Display);
//   glutDisplayFunc(DisplayAll);

    glutMainLoop();

    return 0;
}

void proiectieParalela(unsigned char c) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    switch (c) {
        case 'X':
        case 'Y':
        case 'Z':
            glOrtho(-10, 10, -10, 10, -20, 20);
            break;
        case 'q':
        case 'w':
            glOrtho(-1, 6, -1, 6, -1, 20);
            break;
        default:
            break;
    }
}

void DisplayAxe() {
    int X, Y, Z;
    X = Y = 300;
    Z = 300;
    glLineWidth(2);

    // axa Ox - verde
    glColor3f(0.1, 1, 0.1);
    glBegin(GL_LINE_STRIP);
    glVertex3f(0, 0, 0);
    glVertex3f(X, 0, 0);
    glEnd();

    // axa Oy - albastru
    glColor3f(0.1, 0.1, 1);
    glBegin(GL_LINE_STRIP);
    glVertex3f(0, 0, 0);
    glVertex3f(0, Y, 0);
    glEnd();

    // axa Oz - rosu
    glColor3f(1, 0.1, 0.1);
    glBegin(GL_LINE_STRIP);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, Z);
    glEnd();

    glLineWidth(1);
}

void InitObiect() {
    glNewList(k, GL_COMPILE);
    // fata 1
    glColor3f(1, 0, 0); // rosu
    glBegin(GL_QUADS);
    {
        glVertex3d(0, lat, lat);
        glVertex3d(lat, lat, lat);
        glVertex3d(lat, 0, lat);
        glVertex3d(0, 0, lat);
    }

    glEnd();
    // fata 2
    glColor3f(1, 1, 0); // galben
    glBegin(GL_QUADS);
    {
        glVertex3d(lat, 0, 0);
        glVertex3d(lat, 0, lat);
        glVertex3d(lat, lat, lat);
        glVertex3d(lat, lat, 0);
    }

    glEnd();
    // fata 3
    glColor3f(0, 1, 0); // verde
    glBegin(GL_QUADS);
    {
        glVertex3d(0, lat, lat);
        glVertex3d(lat, lat, lat);
        glVertex3d(lat, lat, 0);
        glVertex3d(0, lat, 0);
    }

    glEnd();
    // fata 4
    glColor3f(0, 0, 1); // albastru
    glBegin(GL_QUADS);
    {
        glVertex3d(0, 0, 0);
        glVertex3d(lat, 0, 0);
        glVertex3d(lat, 0, lat);
        glVertex3d(0, 0, lat);
    }

    glEnd();
    // fata 5
    glColor3f(1, 0, 1);  // magenta
    glBegin(GL_QUADS);
    {
        glVertex3d(0, 0, lat);
        glVertex3d(0, 0, 0);
        glVertex3d(0, lat, 0);
        glVertex3d(0, lat, lat);
    }

    glEnd();
    // fata 6
    glColor3f(0, 1, 1); // cyan
    glBegin(GL_QUADS);
    {
        glVertex3d(0, lat, 0);
        glVertex3d(lat, lat, 0);
        glVertex3d(lat, 0, 0);
        glVertex3d(0, 0, 0);
    }

    glEnd();
    glEndList();
}

void DisplayObiect() {
    glCallList(k);
}