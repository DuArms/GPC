#include <GL/glut.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>

#define or ||
#define and &&

#define x first
#define y second

#define dim 600
#define  PolyDino  (R"(F:\Programming Projects\CLion\Graphic HomeWorks 4\cmake-build-debug\PolyDino.in)")
#define  ExLab  (R"(F:\Programming Projects\CLion\Graphic HomeWorks 4\cmake-build-debug\ex.txt)")

unsigned char prevKey;

const double PI = acos(-1);

template <typename T>
void print(T element){
    std::cout << element << "\n";
}

void drawFilledCircle(double x, double y, double radius, double wRatio, double hRatio){
    int triangleAmount = 20; //# of triangles used to draw circle

    double twicePi = 2.0f * PI;

    glPolygonMode(GL_FRONT, GL_FILL);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y); // center of circle
    for(int i = 0; i <= triangleAmount;i++) {
        glVertex2f(
                x + (radius * std::cos(double(i) *  twicePi / triangleAmount))*wRatio,
                y + (radius * std::sin(double(i) * twicePi / triangleAmount))*hRatio
        );
    }
    glEnd();
}

void drawLine(double x1, double y1, double x2, double y2){
    glVertex2d(x1, y1);
    glVertex2d(x2, y2);
}

struct EdgeBucket {
    int ymax;   //max y-coordinate of edge
    float xofymin;  //x-coordinate of lowest edge point updated only in aet
    float slopeinverse;
};

typedef std::vector<EdgeBucket> EdgeTableTuple;

class FillAlg {
    // start 2
public:
    std::vector<EdgeTableTuple> EdgeTable;
    EdgeTableTuple ActiveEdgeTuple;


// Scanline Function
    void initEdgeTable(int size) {
        EdgeTable.clear();
        for (int i = 0; i < size; i++) {
            EdgeTable.emplace_back(EdgeTableTuple());
        }
    }

    void storeEdgeInTuple(EdgeTableTuple &receiver, int ym, int xm, float slopInv) {

        EdgeBucket tmp{ym, (float) xm, slopInv};


        receiver.push_back(tmp);
        std::sort(ActiveEdgeTuple.begin(), ActiveEdgeTuple.end(),
                  [](EdgeBucket a, EdgeBucket b) -> bool {
                      return a.xofymin < b.xofymin;
                  });

    }

    void storeEdgeInTable(int x1, int y1, int x2, int y2) {
        float m, minv;
        int ymaxTS, xwithyminTS, scanline; //ts stands for to store

        if (x2 == x1) {
            minv = 0.000000;
        } else {
            m = ((float) (y2 - y1)) / ((float) (x2 - x1));

            if (y2 == y1)
                return;

            minv = (float) 1.0 / m;
        }

        if (y1 > y2) {
            scanline = y2;
            ymaxTS = y1;
            xwithyminTS = x2;
        } else {
            scanline = y1;
            ymaxTS = y2;
            xwithyminTS = x1;
        }

        storeEdgeInTuple(EdgeTable[scanline], ymaxTS, xwithyminTS, minv);
    }

    static void removeEdgeByYmax(EdgeTableTuple &Tup, int yy) {
        for (int i = 0; i < Tup.size(); i++) {
            if (Tup[i].ymax == yy) {
                for (int j = i; j < Tup.size() - 1; j++) {
                    Tup[j] = Tup[j + 1];
                }
                Tup.pop_back();
                i--;
            }
        }
    }


    static void updatexbyslopeinv(EdgeTableTuple &Tup) {
        for (auto & i : Tup) {
            i.xofymin = i.xofymin + i.slopeinverse;
        }
    }
};

class CartesianGrid{

private:
    int nrLines{};
    int nrColumns{};

    double yStep{};
    double xStep{};

    double xStart = 0.9;
    double yStart = 0.9;

    double left = -xStart;
    double right = xStart;

    double top = yStart;
    double bot = -yStart;

    double wRatio = 1.0;
    double hRatio = 1.0;

    double wRatioS = 1.0;
    double hRatioS = 1.0;

    double radius{};

    FillAlg fill_alg;

public:
    static double w;
    static double h;

    CartesianGrid(int y, int x) {this->initialize(y, x);}
    explicit CartesianGrid(int x) {this->initialize(x, x);}
    CartesianGrid() {this->initialize(10, 10);}

    inline void resize(double hParam, double wParam){
        w = wParam;
        h = hParam;
        initialize(nrLines, nrColumns);
    }

    inline void initialize(int y, int x){
        nrLines = y;
        nrColumns = x;

        double minRatio = std::min(h/nrLines, w/nrColumns);
        wRatio = minRatio/(w/nrColumns);
        hRatio = minRatio/(h/nrLines);

        left = -xStart * wRatio;
        right = xStart * wRatio;

        top = yStart * hRatio;
        bot = -yStart * hRatio;

        double xSpace = right-left;
        double ySpace = top-bot;

        this->xStep = xSpace / (nrColumns-1);
        this->yStep = ySpace / (nrLines-1);

        radius = std::fmax(xStep, yStep);
        radius /= 2;
        radius *= 0.75;

        minRatio = std::min(h, w);
        wRatioS = minRatio/(w);
        hRatioS = minRatio/(h);

        fill_alg.initEdgeTable(nrLines);
    }

    inline void displayGrid() const{
        glBegin(GL_LINES);
        for(GLdouble x = left; x <= right; x += xStep){
            drawLine(x, bot, x, top);
        }
        drawLine(right, bot, right, top);

        for(GLdouble y = bot; y <= top; y += yStep){
            drawLine(left, y, right, y);
        }
        drawLine(left, top, right, top);
        glEnd();
    }

    void drawPixel(double x, double y, bool swaped=false) const{
        if (swaped){
            std::swap(x, y);
        }

        if (x < 0 or x >= nrColumns)
            return;

        if (y < 0 or y >= nrLines)
            return;


        double py = bot + yStep*y; // y = (py-bot) / yStep
        double px = left + xStep*x; // x = (px-left / xStep

        drawFilledCircle(px, py, radius, wRatioS, hRatioS);
    }

    std::pair<double, double> getCoord(double realX, double realY) const{
        std::pair<double, double> result;

        // [0,w] -> [0,2] -> [-1, 1]
        realX = realX/(w/2)-1;
        realY = realY/(h/2)-1;
        realY *= -1;

        result.y = (realY-bot) / yStep;
        result.x = (realX-left) / xStep;

        return result;
    }

    void drawRealLine(int x1, int y1, int x2, int y2) const{
        if (x1 < 0 or x1 >= nrColumns)
            return;

        if (x2 < 0 or x2 >= nrColumns)
            return;

        if (y1 < 0 or y1 >= nrLines)
            return;

        if (y2 < 0 or y2 >= nrLines)
            return;

        double py1 = bot + yStep*y1;
        double px1 = left + xStep*x1;

        double py2 = bot + yStep*y2;
        double px2 = left + xStep*x2;

        glBegin(GL_LINES);
        drawLine(px1, py1, px2, py2);
        glEnd();
    }

    void drawPixelLine(int x1, int y1, int x2, int y2, int stroke=0) const{
        if (x1 > x2){
            std::swap(x1, x2);
            std::swap(y1, y2);
        }

        bool swaped = false;
        if (abs(y2-y1) > abs(x2-x1)){
            std::swap(x1, y1);
            std::swap(x2, y2);
            swaped = true;
        }

        if (x1 > x2){
            std::swap(x1, x2);
            std::swap(y1, y2);
        }

        int dx = abs(x2-x1);
        int dy = abs(y2-y1);

        int d = 2*dy - dx;
        int dE = 2*dy;
        int dNE = 2*(dy-dx);
        int x = x1;
        int y = y1;

        int sign = 1;
        if (y1 >= y2){
            sign = -1;
        }

        for(int z = -stroke; z <= stroke; z++)
            drawPixel(x, y+z, swaped);
        while (x < x2){
            if (d <= 0){
                d += dE;
                x ++;
            }
            else {
                d += dNE;
                x ++;

                y += sign;
            }
            for(int z = -stroke; z <= stroke; z++)
                drawPixel(x, y+z, swaped);
        }
    }

    std::vector<std::pair<int, int>> getPixelLine(int x1, int y1, int x2, int y2, int stroke=0) const{
        if (x1 > x2){
            std::swap(x1, x2);
            std::swap(y1, y2);
        }

        bool swaped = false;
        if (abs(y2-y1) > abs(x2-x1)){
            std::swap(x1, y1);
            std::swap(x2, y2);
            swaped = true;
        }

        if (x1 > x2){
            std::swap(x1, x2);
            std::swap(y1, y2);
        }

        int dx = abs(x2-x1);
        int dy = abs(y2-y1);

        int d = 2*dy - dx;
        int dE = 2*dy;
        int dNE = 2*(dy-dx);
        int x = x1;
        int y = y1;

        int sign = 1;
        if (y1 >= y2){
            sign = -1;
        }

        std::vector<std::pair<int, int>> result;

        for(int z = -stroke; z <= stroke; z++) {
            if (swaped)
                result.emplace_back(y + z, x);
            else
                result.emplace_back(x, y + z);
            // drawPixel(x, y+z, swaped);
        }
        while (x < x2){
            if (d <= 0){
                d += dE;
                x ++;
            }
            else {
                d += dNE;
                x ++;

                y += sign;
            }
            for(int z = -stroke; z <= stroke; z++) {
                if (swaped)
                    result.emplace_back(y + z, x);
                else
                    result.emplace_back(x, y + z);
                // drawPixel(x, y+z, swaped);
            }
        }
        return result;
    }

    void drawPixelOfCircle(int x, int y, int sx, int sy, bool isElipse=false) const {
        int dx = x-sx;
        int dy = y-sy;

        drawPixel(sx+dx, sy+dy);
        drawPixel(sx-dx, sy-dy);
        drawPixel(sx-dx, sy+dy);
        drawPixel(sx+dx, sy-dy);
        if (x != y and !isElipse){
            drawPixel(sy+dy, sx+dx);
            drawPixel(sy-dy, sx-dx);
            drawPixel(sy-dy, sx+dx);
            drawPixel(sy+dy, sx-dx);
        }
    }

    void drawPixelCircle(int sX, int sY, int R) const{
        int x = sX;
        int y = sY + R;
        int d = 1-R;

        int dE = 3;
        int dSE = -2*R+5;

        drawPixelOfCircle(x, y, sX, sY);
        while (y > x){
            if (d < 0){
                d += dE;
                dE += 2;
                dSE += 2;
            }
            else{
                d += dSE;
                dE += 2;
                dSE += 4;
                y--;
            }
            x++;
            drawPixelOfCircle(x, y, sX, sY);
        }
    }

    void drawPixelCircle2(int sX, int sY, int R, int stroke=0) const{
        int x = sX + R;
        int y = sY;
        int d = 1-R;

        int dE = 3;
        int dSE = -2*R+5;

        for(int z = -stroke; z <= stroke; z++)
            drawPixelOfCircle(x+z, y, sX, sY);
        while (y < x){
            if (d < 0){
                d += dE;
                dE += 2;
                dSE += 2;
            }
            else{
                d += dSE;
                dE += 2;
                dSE += 4;
                x--;
            }
            y++;
            for(int z = -stroke; z <= stroke; z++)
                drawPixelOfCircle(x+z, y, sX, sY);
        }
    }

    void drawPixelElipse(int sX, int sY, int a, int b) const{
        int x = 0;
        int y = b;
        double d1 = b*b - a*a*b + a*a/4.0;
        drawPixelOfCircle(x+sX, y+sY, sX, sY, true);

        while(a*a*(y-0.5) > b*b*(x+1)){
            if (d1 < 0){
                d1 += b*b*(2*x+3);
                x++;
            }
            else{
                d1 += b*b*(2*x+3) + a*a*(-2*y+2);
                x++;
                y--;
            }
            drawPixelOfCircle(x+sX, y+sY, sX, sY, true);
        }

        double d2 = b*b*(x+0.5)*(x+0.5) + a*a*(y-1)*(y-1) - a*a*b*b;

        while (y > 0) {
            if (d2 < 0) {
                d2 +=  b*b*(2*x+2) + a*a*(-2*y+3);
                x++;
                y--;
            }
            else {
                d2 += a*a*(-2*y+3);
                y--;
            }
            drawPixelOfCircle(x+sX, y+sY, sX, sY, true);
        }
    }

    void drawPixelLineEclipse(int y, int x1, int x2, int x0, int y0) const{
        if (x1 > x2)
            std::swap(x1, x2);

        for (int x = x1; x <= x2; x++){
            drawPixelOfCircle(x, y, x0, y0, true);
            // drawPixel(x ,y);
        }
    }

    void drawPixelElipseFill(int x0, int y0, int a, int b) const{
        int xi = 0;
        int x = 0;
        int y = b;
        double fxpyp = 0.0;

        drawPixelLineEclipse(y+y0, x+x0, xi + x0, x0, y0);

        double deltaE, deltaSE, deltaS;

        while (a*a*(y-0.5) > b*b*(x+1)){
            deltaE = b * b * ( 2 * x + 1 );
            deltaSE = b * b * ( 2 * x + 1 ) + a * a * (-2 * y + 1);
            if (fxpyp + deltaE <= 0.0)
            {
                fxpyp += deltaE;
                x++;
                drawPixelLineEclipse(y+y0, x+x0, xi + x0, x0, y0);
            }
            else if (fxpyp + deltaSE <= 0.0)
            {
                fxpyp += deltaSE;
                x++;y--;
                drawPixelLineEclipse(y+y0, x+x0, xi + x0, x0, y0);
            }
        }

        while (y > 0) {
            deltaSE = b * b * (2 * x + 1) + a * a * (-2 * y + 1);
            deltaS = a * a * (-2 * y + 1);
            if (fxpyp + deltaSE <= 0.0) {
                fxpyp += deltaSE;
                x++;
                y--;
            } else {
                fxpyp += deltaS;
                y--;
            }
            drawPixelLineEclipse(y+y0, x+x0, xi + x0, x0, y0);
        }
    }

    void drawPixelElipseFill2(int x0, int y0, int a, int b) const{
        int xi = 0;
        int x = 0;
        int y = -b;
        double fxpyp = 0.0;

        drawPixelLineEclipse(y+y0, x+x0, xi + x0, x0, y0);

        double deltaE, deltaSE, deltaS;

        while (a*a*(y-0.5) < b*b*(x+1)){
            deltaE = b * b * ( 2 * x + 1 );
            deltaSE = b * b * ( 2 * x + 1 ) + a * a * (-2 * y + 1);
            if (fxpyp + deltaE >= 0.0)
            {
                fxpyp += deltaE;
                x--;
                drawPixelLineEclipse(y+y0, x+x0, xi + x0, x0, y0);
            }
            else if (fxpyp + deltaSE >= 0.0)
            {
                fxpyp += deltaSE;
                x--;
                y++;
                drawPixelLineEclipse(y+y0, x+x0, xi + x0, x0, y0);
            }
        }

        while (y < 0) {
            deltaSE = b * b * (2 * x + 1) + a * a * (-2 * y + 1);
            deltaS = a * a * (-2 * y + 1);
            if (fxpyp + deltaSE >= 0.0) {
                fxpyp += deltaSE;
                x--;
                y++;
            } else {
                fxpyp += deltaS;
                y++;
            }
            drawPixelLineEclipse(y+y0, x+x0, xi + x0, x0, y0);
        }
    }

    void drawPolygonFill(const std::vector <std::pair<int, int>>& points) const{
        GLdouble gray = 0.33;

        std::vector <std::pair<int, int>> edgePoints;
        for (int i = 0; i < points.size()-1; i++){
            auto miniResult = getPixelLine(points[i].x, points[i].y, points[i+1].x, points[i+1].y, 0);
            edgePoints.insert(edgePoints.end(), miniResult.begin(), miniResult.end());
            glColor3f(1, 0, 0);
            drawRealLine(points[i].x, points[i].y, points[i+1].x, points[i+1].y);
            glColor3f(gray, gray, gray);
        }
        int n = points.size();
        auto miniResult = getPixelLine(points[n-1].x, points[n-1].y, points[0].x, points[0].y, 0);
        edgePoints.insert(edgePoints.end(), miniResult.begin(), miniResult.end());
        glColor3f(1, 0, 0);
        drawRealLine(points[n-1].x, points[n-1].y, points[0].x, points[0].y);
        glColor3f(gray, gray, gray);


        std::map<std::pair<int, int>, int> edges;
        for (auto const& point : edgePoints){
            drawPixel(point.x, point.y);
            edges[point] = 1;
            //edges[{point.x-1, point.y}] = 1;
            //drawRealLine(13, 13, 10, 0);
        }

        int x_min = nrColumns;
        int y_min = nrLines;
        int x_max = 0;
        int y_max = 0;

        for (auto const& point : points){
            x_min = std::min(x_min, point.x);
            x_max = std::max(x_max, point.x);
            y_min = std::min(y_min, point.y);
            y_max = std::max(y_max, point.y);
        }

        for (int y = y_min; y <= y_max; y++){
            bool outside = true;
            for (int x = x_min-3; x <= x_max+3; x++){
                if (edges[{x, y}] == 1 and edges[{x+1, y}] == 0){
                    outside = !outside;
                }

                if (!outside){
                    //drawPixel(x, y);
                }
            }
        }
    }

    void drawPolyFrame(std::vector<std::pair<int, int>> points) {
        int x1, y1, x2, y2;

        x1 = points[0].x;
        y1 = points[0].y;

        for (int i = 1; i < points.size(); i++) {
            auto p2 = points[i];
            x2 = p2.x;
            y2 = p2.y;
            fill_alg.storeEdgeInTable(x1, y1, x2, y2);//storage of edges in edge table.
            glColor3f(1, 0, 0);
            drawRealLine(x1, y1, x2, y2);
            x1 = x2;
            y1 = y2;
        }
    }

    void drawPolygonFill2(const std::vector<std::pair<int, int>> &points) {
        drawPolyFrame(points);
        ScanlineFill(points, 0, 0.7, 0);
        glFlush();
    }

    void ScanlineFill(const std::vector<std::pair<int, int>> &points, float r, float g, float b) {

        int i, j, x1, ymax1, x2, ymax2, FillFlag = 0, coordCount;

        for (i = 0; i < fill_alg.EdgeTable.size(); i++) {

            for (j = 0; j < fill_alg.EdgeTable[i].size(); j++) {

                fill_alg.storeEdgeInTuple(
                        fill_alg.ActiveEdgeTuple,
                        fill_alg.EdgeTable[i][j].ymax,
                        fill_alg.EdgeTable[i][j].xofymin,
                        fill_alg.EdgeTable[i][j].slopeinverse
                );
            }

            FillAlg::removeEdgeByYmax(fill_alg.ActiveEdgeTuple, i);

            std::sort(fill_alg.ActiveEdgeTuple.begin(), fill_alg.ActiveEdgeTuple.end(),
                      [](EdgeBucket a, EdgeBucket b) -> bool {
                          return a.xofymin < b.xofymin;
                      });

            j = 0;
            FillFlag = 0;
            coordCount = 0;
            x1 = 0;
            x2 = 0;
            ymax1 = 0;
            ymax2 = 0;
            while (j < fill_alg.ActiveEdgeTuple.size()) {
                if (coordCount % 2 == 0) {
                    x1 = (int) (fill_alg.ActiveEdgeTuple[j].xofymin);
                    ymax1 = fill_alg.ActiveEdgeTuple[j].ymax;
                    if (x1 == x2) {
                        if (((x1 == ymax1) && (x2 != ymax2)) || ((x1 != ymax1) && (x2 == ymax2))) {
                            x2 = x1;
                            ymax2 = ymax1;
                        } else {
                            coordCount++;
                        }
                    } else {
                        coordCount++;
                    }
                } else {
                    x2 = (int) fill_alg.ActiveEdgeTuple[j].xofymin;
                    ymax2 = fill_alg.ActiveEdgeTuple[j].ymax;

                    FillFlag = 0;

                    if (x1 == x2) {
                        if (((x1 == ymax1) && (x2 != ymax2)) || ((x1 != ymax1) && (x2 == ymax2))) {
                            x1 = x2;
                            ymax1 = ymax2;
                        } else {
                            coordCount++;
                            FillFlag = 1;
                        }
                    } else {
                        coordCount++;
                        FillFlag = 1;
                    }

                    if (FillFlag) {
                        glColor3f(r, g, b);
                        drawPixelLine(x1, i, x2, i);
                    }
                }
                j++;
            }
            FillAlg::updatexbyslopeinv(fill_alg.ActiveEdgeTuple);
        }
        printf("\nScanline filling complete");
    }
};

double CartesianGrid::w = dim;
double CartesianGrid::h = dim;

CartesianGrid grid;

std::vector<std::pair<double, double>> linesToDraw;

void Init() {
    glClearColor(1.0,1.0,1.0,1.0);
    glLineWidth(1);
    glPointSize(3);
    glPolygonMode(GL_FRONT, GL_LINE);
}

std::vector<std::pair<int, int>> readFromFile(const std::string &name) {
    std::ifstream fin(name);
    std::vector<std::pair<int, int>> buffer;

    int x1, y1;
    while (fin >> x1 >> y1) {
        buffer.emplace_back(x1, y1);
    }

    buffer.push_back(buffer[0]);

    fin.close();
    return buffer;
}

void fillPoly(const std::vector<std::pair<int, int>> &points) {
    //glColor3f(1, gray, gray);
    grid.drawPolygonFill2(points);
}

void Display(){
    auto gray = 0.33f;

    glClear(GL_COLOR_BUFFER_BIT);

    switch (prevKey) {
        case '1': {
            int N = 16;
            grid.initialize(N, 2*N);
            glColor3f(gray, gray, gray);
            grid.displayGrid();

            glColor3f(gray, gray, gray);
            grid.drawPixelLine(0, N-1, N-1, 10, 1);
            grid.drawPixelLine(0, 0, N-1, 7);

            // grid.drawPixelLine(13, 13, 10, 0);

            glColor3f(1, 0, 0);
            grid.drawRealLine(0, N-1, N-1, 10);
            grid.drawRealLine(0, 0, N-1, 7);

            // grid.drawRealLine(13, 13, 10, 0);

//            int i = 0;
//            while (true){
//                if (i+1 >= linesToDraw.size())
//                    break;
//
//                auto p1 = linesToDraw[i];
//                auto p2 = linesToDraw[i+1];
//
//                grid.drawPixelLine(p1.x, p1.y, p2.x, p2.y);
//                grid.drawRealLine(p1.x, p1.y, p2.x, p2.y);
//                i += 2;
//            }
            break;
        }

        case '2': {
            grid.initialize(32, 32);
            glColor3f(gray, gray, gray);
            grid.displayGrid();
            glColor3f(1, gray, gray);
            grid.drawPixelCircle2(16, 16, 12, 1);
            break;
        }

        case '3': {
            grid.initialize(32, 32);
            glColor3f(gray, gray, gray);
            grid.displayGrid();
            glColor3f(1, gray, gray);
            grid.drawPixelElipseFill2(16, 16, 6, 12);
            break;
        }

        case '4': {
            glClear(GL_COLOR_BUFFER_BIT);
            glColor3f(gray, gray, gray);
            grid.initialize(16, 16);
            grid.displayGrid();
            std::vector<std::pair<int, int>> points = readFromFile(ExLab);
            fillPoly(points);
            break;
        }

        case '5': {
            glClear(GL_COLOR_BUFFER_BIT);
            glColor3f(gray, gray, gray);
            grid.initialize(600, 600);
            // grid.displayGrid();
            std::vector<std::pair<int, int>> points = readFromFile(PolyDino);
            fillPoly(points);
            break;
        }

        default:
            break;
    }

    glFlush();
}

void Reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    grid.resize(h, w);
    grid.displayGrid();
    glFlush();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void KeyboardFunc(unsigned char key, int x, int y)
{
    prevKey = key;
    if (key == 27) // escape
        exit(0);
    glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON and state == GLUT_DOWN){
        auto res = grid.getCoord(x, y);
        linesToDraw.push_back(res);
    }

    if (button == GLUT_LEFT_BUTTON and state == GLUT_UP){
        auto res = grid.getCoord(x, y);
        linesToDraw.push_back(res);

        auto p1 = linesToDraw[linesToDraw.size()-2];
        auto p2 = linesToDraw[linesToDraw.size()-1];

        grid.drawPixelLine(p1.x, p1.y, p2.x, p2.y);
        grid.drawRealLine(p1.x, p1.y, p2.x, p2.y);
        glFlush();
    }
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);

    glutInitWindowSize(dim, dim);

    glutInitWindowPosition(100, 100);

    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);

    glutCreateWindow (argv[0]);

    Init();

    glutReshapeFunc(Reshape);

    glutKeyboardFunc(KeyboardFunc);

    glutMouseFunc(MouseFunc);

    glutDisplayFunc(Display);

    glutMainLoop();

    return 0;
}