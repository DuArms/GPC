#include <GL/glut.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
// dimensiunea ferestrei in pixeli
#define dim 300

unsigned char prevKey;

// concoida lui Nicomede (concoida dreptei)
// $x = a + b \cdot cos(t), y = a \cdot tg(t) + b \cdot sin(t)$. sau
// $x = a - b \cdot cos(t), y = a \cdot tg(t) - b \cdot sin(t)$. unde
// $t \in (-\pi / 2, \pi / 2)$
void Display1() {
    double xmax, ymax, xmin, ymin;
    double a = 1, b = 2;
    double pi = 4 * atan(1);
    double ratia = 0.05;
    double t;
    // calculul valorilor maxime/minime ptr. x si y
    // aceste valori vor fi folosite ulterior la scalare
    xmax = a - b - 1;
    xmin = a + b + 1;
    ymax = ymin = 0;
    for (double t = -pi / 2 + ratia; t < pi / 2; t += ratia) {
        double x1, y1, x2, y2;
        x1 = a + b * cos(t);
        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;

        x2 = a - b * cos(t);
        xmax = (xmax < x2) ? x2 : xmax;
        xmin = (xmin > x2) ? x2 : xmin;

        y1 = a * tan(t) + b * sin(t);
        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;

        y2 = a * tan(t) - b * sin(t);
        ymax = (ymax < y2) ? y2 : ymax;
        ymin = (ymin > y2) ? y2 : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    // afisarea punctelor propriu-zise precedata de scalare
    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (t = -pi / 2 + ratia; t < pi / 2; t += ratia) {
        double x1, y1, x2, y2;
        x1 = (a + b * cos(t)) / xmax;
        x2 = (a - b * cos(t)) / xmax;
        y1 = (a * tan(t) + b * sin(t)) / ymax;
        y2 = (a * tan(t) - b * sin(t)) / ymax;

        glVertex2f(x1, y1);
    }
    glEnd();

    glBegin(GL_LINE_STRIP);
    for (t = -pi / 2 + ratia; t < pi / 2; t += ratia) {
        double x1, y1, x2, y2;
        x1 = (a + b * cos(t)) / xmax;
        x2 = (a - b * cos(t)) / xmax;
        y1 = (a * tan(t) + b * sin(t)) / ymax;
        y2 = (a * tan(t) - b * sin(t)) / ymax;

        glVertex2f(x2, y2);
    }
    glEnd();
}

// graficul functiei
// $f(x) = \bar sin(x) \bar \cdot e^{-sin(x)}, x \in \langle 0, 8 \cdot \pi \rangle$,

double functie(double x) {
    return x == 0 ? 1 : abs(round(x) - x) / x;
}

void Display2() {
    double xmax, ymax;
    double ratia = 0.001;
    double t;
    double const min_v = 0 - ratia;
    double const max_v = +100 + ratia;
    xmax = 25;
    ymax = 1.1;

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for (t = min_v + ratia; t < max_v; t += ratia) {
        double x1, y1;
        x1 = t / xmax;
        y1 = functie(t) / ymax;

        glVertex2f(x1, y1);
    }
    glEnd();
}


#include <vector>

std::vector<double>
get_function_scale(double d_min, double d_max, double rate, double (*function_x)(double),
                   double (*function_y)(double)) {
    double x1, y1;
    double xmax, ymax, xmin, ymin;

    xmax = xmin = function_x(d_min + rate);
    ymax = ymin = function_y(d_min + rate);

    for (double t = d_min + 2 * rate; t < d_max; t += rate) {
        x1 = function_x(t);
        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;

        y1 = function_y(t);
        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    std::vector<double> rez;
    rez.push_back(xmax);
    rez.push_back(ymax);

    return rez;
}

double x_pascal(double t) { return (2 * (0.3 * cos(t) + 0.2) * cos(t)); }

double y_pascal(double t) { return (2 * (0.3 * cos(t) + 0.2) * sin(t)); }

void Display3() {
    double xmax, ymax;
    double const pi = acos(-1);
    double ratia = 0.05;
    double t;

    auto values = get_function_scale(-pi + ratia, pi, ratia, x_pascal, y_pascal);
    xmax = values[0] + 0.2;
    ymax = values[1] + 0.2;

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for (t = -pi + ratia; t < pi; t += ratia) {
        double x1, y1;
        x1 = x_pascal(t) / xmax;
        y1 = y_pascal(t) / ymax;

        glVertex2f(x1, y1);
    }
    glVertex2f(x_pascal(-pi + ratia) / xmax, y_pascal(-pi + ratia) / ymax);

    glEnd();
}

double x_trisectoarea_lui_Longchamps(double t) { return 0.2 / (4 * cos(t) * cos(t) - 3); }

double y_trisectoarea_lui_Longchamps(double t) { return 0.2 * tan(t) / (4 * cos(t) * cos(t) - 3); }

void Display4() {
    double xmax, ymax;
    double const pi = acos(-1);
    double ratia = 0.025;
    double t;
    auto values = get_function_scale(-pi / 2, pi / 2, ratia, x_trisectoarea_lui_Longchamps,
                                     y_trisectoarea_lui_Longchamps);

    xmax = values[0] / 3.5;
    ymax = values[1] / 3.5;

    std::cout << xmax << " " << ymax;

    glColor3f(0.1, 0.1, 1); // rosu
    glBegin(GL_LINE_STRIP);
#define x_scale (2)
    double x1, y1, x2, y2;
    double x_colt, y_colt;
    int h((-pi / 6 - ratia + pi / 2) / ratia);
    double l(-pi / 2 + ratia * (h + 1));

    x1 = (x_trisectoarea_lui_Longchamps(l) / xmax) * x_scale;
    y1 = y_trisectoarea_lui_Longchamps(l) / ymax;

    x2 = (x_trisectoarea_lui_Longchamps(-pi / 2 + ratia) / xmax) * x_scale;
    y2 = y_trisectoarea_lui_Longchamps(-pi / 2 + ratia) / ymax;

    x_colt = x1;
    y_colt = y2;

    glVertex2f(x1, y1);
    glVertex2f(x1, y2);
    glVertex2f(x2, y2);

    std::vector<std::pair<double, double>> poz_triunghi;

    for (t = -pi / 2 + ratia; t < -pi / 6; t += ratia) {
        if (abs(t) == pi / 6) {
            continue;
        }

        x1 = (x_trisectoarea_lui_Longchamps(t) / xmax) * x_scale;
        y1 = y_trisectoarea_lui_Longchamps(t) / ymax;

        for (double r = t; r < t + ratia && r < -pi / 6 - ratia; r += ratia / 8) {
            poz_triunghi.emplace_back((x_trisectoarea_lui_Longchamps(r) / xmax) * x_scale,
                                      y_trisectoarea_lui_Longchamps(r) / ymax);
        }
        glVertex2f(x1, y1);
    }

    glEnd();
    glColor3f(1, 0.1, 0.1); // rosu
    glPolygonMode(GL_FRONT, GL_FILL);
    glBegin(GL_TRIANGLES);

    for (int i = 0; i < poz_triunghi.size() - 1; i += 3) {
        x1 = poz_triunghi[i].first;
        y1 = poz_triunghi[i].second;

        x2 = poz_triunghi[i + 1].first;
        y2 = poz_triunghi[i + 1].second;

        glVertex2f(x1, y1);
        glVertex2f(x_colt, y_colt);
        glVertex2f(x2, y2);
    }
    glEnd();
}

double x_cicloida(double t) { return 0.1 * t - 0.2 * sin(t); }

double y_cicloida(double t) { return 0.1 - 0.2 * cos(t); }

void Display5() {
    double xmax, ymax;
    double const pi = acos(-1);
    double ratia = 0.1;
    double t;
    double const min_v = -1000;
    double const max_v = +1000;

    auto values = get_function_scale(min_v, max_v, ratia, x_cicloida,
                                     y_cicloida);
    xmax = 1;
    ymax = 1;

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for (t = min_v + ratia; t < max_v; t += ratia) {
        double x1, y1;
        x1 = x_cicloida(t) / xmax;
        y1 = y_cicloida(t) / ymax;

        glVertex2f(x1, y1);
    }
    glEnd();

}

double x_epicicloida(double t) {
    const double R = 0.1, r = 0.3;
    return (R + r) * cos(r / R * t) - r * cos(t + r / R * t);
}

double y_epicicloida(double t) {
    const double R = 0.1, r = 0.3;
    return (R + r) * sin(r / R * t) - r * sin(t + r / R * t);
}

void Display6() {
    double xmax, ymax;
    double const pi = acos(-1);
    double ratia = 0.05;
    double t;
    double const min_v = 0;
    double const max_v = 2 * pi;

    auto values = get_function_scale(min_v, max_v, ratia, x_epicicloida,
                                     y_epicicloida);
    xmax = values[0] + 0.2;
    ymax = values[1] + 0.2;

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for (t = min_v + ratia; t < max_v; t += ratia) {
        double x1, y1;
        x1 = x_epicicloida(t) / xmax;
        y1 = y_epicicloida(t) / ymax;

        glVertex2f(x1, y1);
    }
    glEnd();
}

double x_hipocicloida(double t) {
    const double R = 0.1, r = 0.3;
    return (R - r) * cos(r / R * t) - r * cos(t - r / R * t);
}

double y_hipocicloida(double t) {
    const double R = 0.1, r = 0.3;
    return (R - r) * sin(r / R * t) - r * sin(t - r / R * t);
}

void Display7() {
    double xmax, ymax;
    double const pi = acos(-1);
    double ratia = pi / 180;
    double t;
    double const min_v = 0;
    double const max_v = 2 * pi;

    auto values = get_function_scale(min_v, max_v, ratia, x_hipocicloida,
                                     y_hipocicloida);
    xmax = values[0] + 0.2;
    ymax = values[1] + 0.2;

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for (t = min_v + ratia; t < max_v; t += ratia) {
        double x1, y1;
        x1 = x_hipocicloida(t) / xmax;
        y1 = y_hipocicloida(t) / ymax;

        glVertex2f(x1, y1);
    }
    glEnd();
}

double x_polar(double r, double t) {
    return r * cos(t);
}

double y_polar(double r, double t) {
    return r * sin(t);
}

double r_lemniscata_lui_Bernoulli(double t) {
    return 0.4 * sqrt(2 * cos(2 * t));
}

void Display8() {
    double xmax, ymax;
    double const pi = acos(-1);
    double ratia = pi / 480;
    double t;
    double const min_v = -pi / 4;
    double const max_v = pi / 4;

//    auto values = get_function_scale(min_v, max_v, ratia, x_hipocicloida,
//                                     y_hipocicloida);
    xmax = 0.75;
    ymax = 0.75;

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for (t = max_v - ratia; t > min_v; t -= ratia) {
        double x1, y1;
        double r = r_lemniscata_lui_Bernoulli(t);

        x1 = x_polar(r, t) / xmax;
        y1 = y_polar(r, t) / ymax;

        glVertex2f(x1, y1);


    }
    for (t = min_v + ratia; t < max_v; t += ratia) {
        double x1, y1;
        double r = -r_lemniscata_lui_Bernoulli(t);

        x1 = x_polar(r, t) / xmax;
        y1 = y_polar(r, t) / ymax;

        glVertex2f(x1, y1);
    }
    glEnd();
}

double r_spirala_logaritmica(double t) {
    return 0.02 * exp(t + 1);
}

void Display9() {
    double xmax, ymax;
    double const pi = acos(-1);
    double ratia = 0.05;
    double t;
    double const min_v = 0;
    double const max_v = 1000;

    xmax = 50;
    ymax = 50;

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for (t = min_v + ratia; t < max_v; t += ratia) {
        double x1, y1;
        double r = -r_spirala_logaritmica(t);

        x1 = x_polar(r, t) / xmax;
        y1 = y_polar(r, t) / ymax;

        glVertex2f(x1, y1);
    }
    glEnd();
}

void Init(void) {

    glClearColor(1.0, 1.0, 1.0, 1.0);

    glLineWidth(1);

//   glPointSize(4);

    glPolygonMode(GL_FRONT, GL_LINE);
}

void Display(void) {
    glClear(GL_COLOR_BUFFER_BIT);

    switch (prevKey) {
        case '1':
            Display1();
            break;
        case '2':
            Display2();
            break;
        case '3':
            Display3();
            break;
        case '4':
            Display4();
            break;
        case '5':
            Display5();
            break;
        case '6':
            Display6();
            break;
        case '7':
            Display7();
            break;
        case '8':
            Display8();
            break;
        case '9':
            Display9();
            break;

        default:
            break;
    }

    glFlush();
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