#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include <Gl/glut.h>

// dimensiunea ferestrei in pixeli
#define dim 600
// numarul maxim de iteratii pentru testarea apartenentei la mult.Julia-Fatou
#define NRITER_JF 5000
// modulul maxim pentru testarea apartenentei la mult.Julia-Fatou
#define MODMAX_JF 10000000
// ratii ptr. CJuliaFatou
#define RX_JF 0.01
#define RY_JF 0.01

unsigned char prevKey;

class CComplex {
public:
    CComplex() : re(0.0), im(0.0) {}

    CComplex(double re1, double im1) : re(re1 * 1.0), im(im1 * 1.0) {}

    CComplex(const CComplex &c) : re(c.re), im(c.im) {}

    ~CComplex() {}

    CComplex &operator=(const CComplex &c) {
        re = c.re;
        im = c.im;
        return *this;
    }

    double getRe() { return re; }

    void setRe(double re1) { re = re1; }

    double getIm() { return im; }

    void setIm(double im1) { im = im1; }

    double getModul() { return sqrt(re * re + im * im); }

    int operator==(CComplex &c1) {
        return ((re == c1.re) && (im == c1.im));
    }


    CComplex pow2() {
        CComplex rez;
        rez.re = powl(re * 1.0, 2) - powl(im * 1.0, 2);
        rez.im = 2.0 * re * im;
        return rez;
    }

    friend CComplex operator+(const CComplex &c1, const CComplex &c2);

    friend CComplex operator*(CComplex &c1, CComplex &c2);

    bool operator<(CComplex &c2) {
        return this->getModul() < c2.getModul();
    }

    void print(FILE *f) {
        fprintf(f, "%.20f%+.20f i", re, im);
    }

private:
    double re, im;
};

CComplex operator+(const CComplex &c1, const CComplex &c2) {
    CComplex rez(c1.re + c2.re, c1.im + c2.im);
    return rez;
}

CComplex operator*(CComplex &c1, CComplex &c2) {
    CComplex rez(c1.re * c2.re - c1.im * c2.im,
                 c1.re * c2.im + c1.im * c2.re);
    return rez;
}

class CJuliaFatou {
public:
    CJuliaFatou() {
        // m.c se initializeaza implicit cu 0+0i

        m.nriter = NRITER_JF;
        m.modmax = MODMAX_JF;
    }

    CJuliaFatou(CComplex &c) {
        m.c = c;
        m.nriter = NRITER_JF;
        m.modmax = MODMAX_JF;
    }

    ~CJuliaFatou() {}

    void setmodmax(double v) {
        assert(v <= MODMAX_JF);
        m.modmax = v;
    }

    double getmodmax() { return m.modmax; }

    void setnriter(int v) {
        assert(v <= NRITER_JF);
        m.nriter = v;
    }

    int getnriter() { return m.nriter; }

    // testeaza daca x apartine multimii Julia-Fatou Jc
    // returneaza 0 daca apartine, -1 daca converge finit, +1 daca converge infinit
    int isIn(CComplex &x) {
        int rez = 0;
        // tablou in care vor fi memorate valorile procesului iterativ z_n+1 = z_n * z_n + c
        CComplex z0, z1;

        z0 = x;
        for (int i = 1; i < m.nriter; i++) {
            z1 = z0 * z0 + m.c;
            if (z1 == z0) {
                // x nu apartine m.J-F deoarece procesul iterativ converge finit
                rez = -1;
                break;
            } else if (z1.getModul() > m.modmax) {
                // x nu apartine m.J-F deoarece procesul iterativ converge la infinit
                rez = 1;
                break;
            }
            z0 = z1;
        }

        return rez;
    }

    // afisarea multimii J-F care intersecteaza multimea argument
    void display(double xmin, double ymin, double xmax, double ymax) {
        glPushMatrix();
        glLoadIdentity();

//    glTranslated((xmin + xmax) * 1.0 / (xmin - xmax), (ymin + ymax)  * 1.0 / (ymin - ymax), 0);
//    glScaled(1.0 / (xmax - xmin), 1.0 / (ymax - ymin), 1);
        // afisarea propriu-zisa
        glBegin(GL_POINTS);
        for (double x = xmin; x <= xmax; x += RX_JF)
            for (double y = ymin; y <= ymax; y += RY_JF) {
                CComplex z(x, y);
                int r = isIn(z);
//        z.print(stdout);
                if (r == 0) {
//          fprintf(stdout, "   \n");
                    glVertex3d(x, y, 0);
                } else if (r == -1) {
//          fprintf(stdout, "   converge finit\n");
                } else if (r == 1) {
//          fprintf(stdout, "   converge infinit\n");
                }
            }
        fprintf(stdout, "STOP\n");
        glEnd();

        glPopMatrix();
    }

private:
    struct SDate {
        CComplex c;
        // nr. de iteratii
        int nriter;
        // modulul maxim
        double modmax;
    } m;
};

// multimea Julia-Fatou pentru z0 = 0 si c = -0.12375+0.056805i
void Display1() {
    CComplex c(-0.12375, 0.056805);
    CJuliaFatou cjf(c);

    glColor3f(1.0, 0.1, 0.1);
    cjf.setnriter(30);
    cjf.display(-0.8, -0.4, 0.8, 0.4);
}

// multimea Julia-Fatou pentru z0 = 0 si c = -0.012+0.74i
void Display2() {
    CComplex c(-0.012, 0.74);
    CJuliaFatou cjf(c);

    glColor3f(1.0, 0.1, 0.1);
    cjf.setnriter(30);
    cjf.display(-1, -1, 1, 1);
}
// cod V

#include <iostream>
#include <ctype.h>
#include <vector>
#include <set>
#include <map>


#define my_point_size (1.5)
#define my_step (0.007)


class Multimea_Mandelbrot {
public:
    static const unsigned int n_max = 100;
    static const unsigned long infinite = 2;

    Multimea_Mandelbrot() = default;

    static int check_if_valid_number(const CComplex& complex_number) {

        CComplex z_curent(0, 0);

        for (uint32_t n = 1; n < n_max; n++) {
            z_curent = z_curent.pow2() + complex_number;

            if (z_curent.getModul() >= infinite) {
                return n;
            }
        }
        return n_max;
    }
};

void Display3() {
    glPushMatrix();
    glLoadIdentity();
    glPointSize(my_point_size);
    const double pas = my_step;
    CComplex start_point(0, 0);
    Multimea_Mandelbrot mm;

    glBegin(GL_POINTS);
    for (double x = -2; x < 1; x += pas)
        for (double y = -1; y < 1; y += pas) {
            start_point.setRe(x);
            start_point.setIm(y);
            const int steps = mm.check_if_valid_number(start_point);
            if (steps != Multimea_Mandelbrot::n_max) {
                double ration = (double) steps / Multimea_Mandelbrot::n_max;

                long culoare = (long) (ration * 0xFFFFFF);
                short  blue  = (culoare) & 0xFF;
                short red = (culoare >> 8) & 0xFF;
                short green = (culoare >> 16) & 0xFF;

                glColor3f(red / 255.0f, green / 255.0f, blue / 255.0f);
                glVertex3d((x + 2) * 0.6666 - 1, y * 0.90, 0);
            }
        }
    glEnd();
}

void Display4() {
    glPushMatrix();
    glLoadIdentity();
    glPointSize(my_point_size);
    const float pas = my_step;

    CComplex start_point(0, 0);
    Multimea_Mandelbrot mm;

    glColor3f(1.0, 0.1, 0.1);

    glBegin(GL_POINTS);

    for (double x = -2; x < 1; x += pas)
        for (double y = -1; y < 1; y += pas) {
            start_point.setRe(x);
            start_point.setIm(y);
            if (mm.check_if_valid_number(start_point) == Multimea_Mandelbrot::n_max) {
                glVertex3d((x + 2) * 0.6666 - 1, y * 0.90, 0);
            }
        }
    glEnd();
}

// ex 4

#define dim 300

int nivel = 0;

class C2coord {
public:
    C2coord() {
        m.x = m.y = 0;
    }

    C2coord(double x, double y) {
        m.x = x;
        m.y = y;
    }

    C2coord(C2coord &p) {
        m.x = p.m.x;
        m.y = p.m.y;
    }

    C2coord &operator=(C2coord &p) {
        m.x = p.m.x;
        m.y = p.m.y;
        return *this;
    }

    int operator==(C2coord &p) {
        return ((m.x == p.m.x) && (m.y == p.m.y));
    }

protected:
    struct SDate {
        double x, y;
    } m;
};

class CPunct : public C2coord {
public:
    CPunct() : C2coord(0.0, 0.0) {}

    CPunct(double x, double y) : C2coord(x, y) {}

    CPunct &operator=(const CPunct &p) {
        m.x = p.m.x;
        m.y = p.m.y;
        return *this;
    }

    void getxy(double &x, double &y) {
        x = m.x;
        y = m.y;
    }

    int operator==(CPunct &p) {
        return ((m.x == p.m.x) && (m.y == p.m.y));
    }

    void marcheaza() {
        glBegin(GL_POINTS);
        glVertex2d(m.x, m.y);
        glEnd();
    }

    void print(FILE *fis) {
        fprintf(fis, "(%+f,%+f)", m.x, m.y);
    }

};

double r, g, b;

class CVector : public C2coord {
public:
    CVector() : C2coord(0.0, 0.0) {
        normalizare();
    }

    CVector(double x, double y) : C2coord(x, y) {
        normalizare();
    }

    CVector &operator=(CVector &p) {
        m.x = p.m.x;
        m.y = p.m.y;
        return *this;
    }

    int operator==(CVector &p) {
        return ((m.x == p.m.x) && (m.y == p.m.y));
    }

    double getAngle() {
        return atan2(m.y, m.x);
    }

    CPunct getDest(CPunct &orig, double lungime) {
        double x, y;
        orig.getxy(x, y);
        CPunct p(x + m.x * lungime, y + m.y * lungime);
        return p;
    }

    void rotatie(double grade) {
        double x = m.x;
        double y = m.y;
        double t = 2 * (4.0 * atan(1)) * grade / 360.0;
        m.x = x * cos(t) - y * sin(t);
        m.y = x * sin(t) + y * cos(t);
        normalizare();
    }

    void full_r(double grade) {
        double x = 1;
        double y = 0;

        double t = 2 * (4.0 * atan(1)) * grade / 360.0;
        m.x = x * cos(t) - y * sin(t);
        m.y = x * sin(t) + y * cos(t);

        normalizare();
    }

    void deseneaza(CPunct p, double lungime) {
        double x, y;
        p.getxy(x, y);
        glColor3f(r, g, b);
        glBegin(GL_LINE_STRIP);
//
//        x = x * 2 + 1;
//        y = y * 2 - 1;
//        lungime *= 2;

        glVertex2d(x, y);
        glVertex2d(x + m.x * lungime, y + m.y * lungime);
        glEnd();
    }

    void print(FILE *fis) {
        fprintf(fis, "%+fi %+fj", C2coord::m.x, C2coord::m.y);
    }

private:
    void normalizare() {
        double d = sqrt(C2coord::m.x * C2coord::m.x + C2coord::m.y * C2coord::m.y);
        if (d != 0.0) {
            C2coord::m.x = C2coord::m.x * 1.0 / d;
            C2coord::m.y = C2coord::m.y * 1.0 / d;
        }
    }
};

class MengerSponge {
public:
    void deseneaza_menger_sponge(double length, int nivel, CPunct &p, CVector &v, int d) {
        if (nivel == 0) return;

        double x, y;

        p.getxy(x, y);

        int coord_x[] = {-1, 0, 1, -1, 1, -1, 0, 1};
        int coord_y[] = {-1, -1, -1, 0, 0, 1, 1, 1};
        double new_length = length;

        for (int i = 0; i < 8; i++) {
            CPunct new_p(x + coord_x[i] * new_length, y + coord_y[i] * new_length);
            deseneaza_menger_sponge(new_length / 3, nivel - 1, new_p, v, d);

        }

        //desenez patratul
        p = CPunct(x + length / 2, y + length / 2);

        for (int i = 0; i < 4; i++) {
            v.rotatie(d * 90);
            v.deseneaza(p, length);

            p = v.getDest(p, length);
        }


    }


    void afisare(double lungime, int nivel) {
        CVector v(0.0, 1);
        CPunct p(0, 0);

        r = 1;
        g = 0.1;
        b = 0.1;
        deseneaza_menger_sponge(lungime, nivel, p, v, 1);

        lungime *= 3;
        p = CPunct(lungime / 2, lungime / 2);
        for (int i = 0; i < 4; i++) {
            v.rotatie(1 * 90);
            v.deseneaza(p, lungime);

            p = v.getDest(p, lungime);
        }
    }
};

void Display5() {
    MengerSponge img1;
    img1.afisare(0.5, nivel + 2);
    std::cout << nivel << std::endl;
    nivel++;
}

class CopacCuProbleme {
public:
    void copac_cu_probleme(double length, int nivel, CPunct &p, CVector &v, int d) {
        if (nivel <= 0) return;

        CPunct original_point = p;
        CPunct new_p;
        CVector new_v;

        v.deseneaza(p, length / 4);
        p = v.getDest(p, length / 4);

        v.rotatie(-45);
        v.deseneaza(p, length);

        aux_recursion(length, nivel, d, p, v, new_p, new_v);

        v.rotatie(90);
        v.deseneaza(p, length);

        p = v.getDest(p, length);

        v.rotatie(15);
        v.deseneaza(p, length);

        aux_recursion(length, nivel, d, p, v, new_p, new_v);

        v.rotatie(-15 - 45);
        v.deseneaza(p, length);

        p = v.getDest(p, length);

        v.rotatie(+30);
        v.deseneaza(p, length / 2);

        aux_recursion(length / 2, nivel, d, p, v, new_p, new_v);

        v.rotatie(-120);
        v.deseneaza(p, length / 2);

        aux_recursion(length / 2, nivel, d, p, v, new_p, new_v);
    }

    void aux_recursion(double length, int nivel, int d, CPunct &p, CVector &v, CPunct &new_p, CVector &new_v) {
        new_p = v.getDest(p, length );
        new_v = v;

        std::swap(r, b);
        copac_cu_probleme(length / 2, nivel - 1, new_p, new_v, d);
        std::swap(r, b);
    }

    void afisare(double lungime, int nivel) {
        CVector v(0.0, -1);
        CPunct p(-0.2, 0.80);
        r = 1;
        g = 0.1;
        b = 0.1;
        copac_cu_probleme(0.33, nivel, p, v, 1);
    }
};

void Display6() {
    CopacCuProbleme copac;
    copac.afisare(0.5, nivel);
    std::cout << nivel;
    nivel++;
}

class SierpinskiTriangle {
public:
    void sierpinski_triangle(double length, int nivel, CPunct &p, CVector &v, int d) {
        if (nivel <= 0) return;
        if (nivel == 1) {

            CVector new_v = v;
            int teta;
            int omega;

            if (d == 1) {
                teta = 60;
                omega = -30;
            } else {
                teta = -60;
                omega = +30 + 180;
            }
            v.rotatie(omega);

            for (int i = 0; i < 3; i++) {
                v.rotatie(teta);
                v.deseneaza(p, length);
                p = v.getDest(p, length);
            }
            v = new_v;

            return;
        }
        CVector old_v = v;
        v.rotatie(-d * 60);
        sierpinski_triangle(length * 0.55, nivel - 1, p, v, -d);


        v.rotatie(d * 60);
        sierpinski_triangle(length * 0.55, nivel - 1, p, v, d);

        v.rotatie(d * 60);

        sierpinski_triangle(length * 0.55, nivel - 1, p, v, -d);

        v = old_v;
    }

    void afisare(double lungime, int nivel) {
        CVector v(1.0, 0.0);
        CPunct p(-0.8, -0.8);

        r = 1;
        g = 0.1;
        b = 0.1;

        sierpinski_triangle(lungime, nivel, p, v, 1);
    }

};


void Display7() {
    SierpinskiTriangle triunghi;
    triunghi.afisare(0.35   , nivel);
    std::cout << nivel;
    nivel++;
}

void Display(void) {
    switch (prevKey) {
        case '0':
            glClear(GL_COLOR_BUFFER_BIT);
            nivel = 0;
            fprintf(stderr, "nivel = %d\n", nivel);
            break;
        case '1':
            glClear(GL_COLOR_BUFFER_BIT);
            Display1();
            break;
        case '2':
            glClear(GL_COLOR_BUFFER_BIT);
            Display2();
            break;
        case '3':
            glClear(GL_COLOR_BUFFER_BIT);
            Display3();
            break;
        case '4':
            glClear(GL_COLOR_BUFFER_BIT);
            Display4();
            break;
        case '5':
            glClear(GL_COLOR_BUFFER_BIT);
            Display5();
            break;
        case '6':
            glClear(GL_COLOR_BUFFER_BIT);
            Display6();
            break;
        case '7':
            glClear(GL_COLOR_BUFFER_BIT);
            Display7();
            break;
        default:
            break;
    }

    glFlush();
}













void Init(void) {
    r = 1;
    g = 0.1;
    b = 0.1;
    glClearColor(1.0, 1.0, 1.0, 1.0);

    glLineWidth(1);

//   glPointSize(3);

    glPolygonMode(GL_FRONT, GL_LINE);
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

int main(int argc, char **argv) {

    glutInit(&argc, argv);

    glutInitWindowSize(dim, dim);

    glutInitWindowPosition(100, 100);

    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);

    glutCreateWindow(argv[0]);

    Init();

    glutReshapeFunc(Reshape);

    glutKeyboardFunc(KeyboardFunc);

    glutMouseFunc(MouseFunc);

    glutDisplayFunc(Display);

    glutMainLoop();

    return 0;
}


