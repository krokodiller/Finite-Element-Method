#pragma once
#include <vector>
#include "FiniteElement.h"
class ForwardProblem
{
public:

    ForwardProblem();
    inline double u(double x, double y) { return x * x * x * y * y * y; }
    inline double f(double x, double y, double lambda, double gamma) { return -6 * x * y * y * y * lambda + -6 * y * x * x * x * lambda + gamma * u(x, y); }
    //inline double u(double x, double y) { return x * x * x * y * y; }
    //inline double f(double x, double y, double lambda, double gamma) { return -6 * x * y * y * lambda + -2 * x * x * lambda + gamma * u(x, y); }
    //inline double u(double x, double y) { return 0; }
    //inline double f(double x, double y, double lambda, double gamma) { return 0 * lambda + gamma * u(x, y); }

    void setGridx(std::vector<double> grid);
    void setGridx(double x0, double x1, double step);
    void setGridx(double x0, double x1, int numberOfFes, double coef);
    void setGridy(std::vector<double> grid);
    void setGridy(double y0, double y1, double step);
    void setGridy(double y0, double y1, int numberOfFes, double coef);
    void addArea(double xleftbottom, double yleftbottom, double xrighttop, double yrighttop, double lambda, double gamma);
    void addArea(double xleftbottom, double yleftbottom, double xrighttop, double yrighttop, double sigma);
    void setSourcesPower(double power);
    void setSourcesLocation(double x1, double y1, double x2, double y2);

    void assembleGrid();

    void clear();

    std::vector<double> getMainMatrix();

    std::vector<double> getGMatrix();
    std::vector<double> getMMatrix();
    std::vector<double> getBVector();

    std::vector<double> getQVector();

    std::vector<double> getTrueSolution();

    double getSolutionInPoint(double x, double y);

    double getDerivativeInPoint(double x, double y);

    Point getGrad(double x, double y);

    void solveProblem();


private:

    void applyFirstBorder();
    void applySources();

    double bifunc2d(int number, Point p, int currentFE);
    double bifuncx(int number, double x, int currentFE);
    double bifuncy(int number, double y, int currentFE);

    Point gradbifunc2d(int number, Point p, int currentFE);
    double dbifuncx(int number, double x, int currentFE);
    double dbifuncy(int number, double y, int currentFE);

    void calculateGMatrix();
    void calculateMMatrix();
    void calculateBVector();
    void calculateQVector();

    std::vector<double> getLocalGMatrix(int currentFE);
    std::vector<double> getLocalCMatrix(int currentFE);
    std::vector<double> getLocalMMatrix(int currentFE);
    std::vector<double> getLocalBVector(int currentFE);



    std::vector<int> getGlobalPositionMatrix(int currentFE);
    std::vector<int> getGlobalPositionVector(int currentFE);

    double getLambda(int currentFE);
    double getGamma(int currentFE);
    double getSigma(int currentFE);

    double power = 0;
    Point source0;
    Point source1;

    std::vector<Area> areas;

    std::vector<FiniteElement> fes;
    std::vector<Node> nodes;

    std::vector<double> gridx, gridy;

    std::vector<double> mainMatrix, GMatrix, MMatrix, bVector, qVector;

    int numberOfFiniteElementsx = 0, numberOfFiniteElementsy = 0;
    int globalMatrixRank = 0;
};

