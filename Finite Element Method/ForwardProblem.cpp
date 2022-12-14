#include "ForwardProblem.h"
#include "Solver.h"
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <limits>
constexpr double eps = 0.00001;

ForwardProblem::ForwardProblem()
{

}

void ForwardProblem::setGridx(std::vector<double> grid)
{
    if (grid.empty())
    {
        throw std::invalid_argument("grid is empty");
    }
    if (!std::is_sorted(grid.begin(), grid.end()))
    {
        std::cerr << "grid is not sorted. Sorting...\n";
        std::sort(grid.begin(), grid.end());
    }
    if (std::adjacent_find(grid.begin(), grid.end()) != grid.end())
    {
        std::cerr << "grid values are not uniqe. Removing duplicates...\n";
        auto last = std::unique(grid.begin(), grid.end());
        grid.erase(last, grid.end());
    }
    this->gridx = grid;
    numberOfFiniteElementsx = grid.size() - 1;
}

void ForwardProblem::setGridy(std::vector<double> grid)
{
    if (grid.empty())
    {
        throw std::invalid_argument("grid is empty");
    }
    if (!std::is_sorted(grid.begin(), grid.end()))
    {
        std::cerr << "grid is not sorted. Sorting...\n";
        std::sort(grid.begin(), grid.end());
    }
    if (std::adjacent_find(grid.begin(), grid.end()) != grid.end())
    {
        std::cerr << "grid values are not uniqe. Removing duplicates...\n";
        auto last = std::unique(grid.begin(), grid.end());
        grid.erase(last, grid.end());
    }
    this->gridy = grid;
    numberOfFiniteElementsy = grid.size() - 1;
}

void ForwardProblem::setGridx(double x0, double x1, double step)
{
    gridx.clear();
    int i = 0;
    for (double x = x0; x <= x1 + eps; x = x0 + step * i)
    {
        gridx.push_back(x);
        i++;
    }
    numberOfFiniteElementsx = gridx.size() - 1;
}

void ForwardProblem::setGridx(double x0, double x1, int numberOfFes, double coef)
{
    if (coef < 0)
    {
        coef = 1 / coef;
    }
    double length = x1 - x0;
    double h;
    if (fabs(fabs(coef) - 1) < eps)
    {
        h = length / numberOfFes;
    }
    else
    {
        h = length * ((1 - coef) / (1 - pow(coef, numberOfFes)));
    }


    gridx.push_back(x0);
    double newx = x0;
    for (size_t i = 0; i < numberOfFes; i++)
    {
        newx += (x1 - x0) / length * h;
        h = h * coef;
        gridx.push_back(newx);
    }
    numberOfFiniteElementsx = gridx.size() - 1;
}

void ForwardProblem::setGridy(double y0, double y1, double step)
{
    gridy.clear();
    int i = 0;
    for (double y = y0; y <= y1 + eps; y = y0 + step * i)
    {
        gridy.push_back(y);
        i++;
    }
    numberOfFiniteElementsy = gridy.size() - 1;
}

void ForwardProblem::setGridy(double y0, double y1, int numberOfFes, double coef)
{
    if (coef < 0)
    {
        coef = 1 / coef;
    }
    double length = y1 - y0;
    double h;
    if (fabs(fabs(coef) - 1) < eps)
    {
        h = length / numberOfFes;
    }
    else
    {
        h = length * ((1 - coef) / (1 - pow(coef, numberOfFes)));
    }


    gridy.push_back(y0);
    double newy = y0;
    for (size_t i = 0; i < numberOfFes; i++)
    {
        newy += (y1 - y0) / length * h;
        h = h * coef;
        gridy.push_back(newy);
    }
    numberOfFiniteElementsy = gridy.size() - 1;
}

void ForwardProblem::addArea(double xleftbottom, double yleftbottom, double xrighttop, double yrighttop, double lambda, double gamma)
{
    areas.push_back(Area(xleftbottom, yleftbottom, xrighttop, yrighttop, lambda, gamma));
}

void ForwardProblem::addArea(double xleftbottom, double yleftbottom, double xrighttop, double yrighttop, double sigma)
{
    areas.push_back(Area(xleftbottom, yleftbottom, xrighttop, yrighttop, sigma));
}

void ForwardProblem::setSourcesPower(double power)
{
    this->power = power;
}

void ForwardProblem::setSourcesLocation(double x1, double y1, double x2, double y2)
{
    source0.x = x1;
    source0.y = y1;
    source1.x = x2;
    source1.y = y2;
}

void ForwardProblem::assembleGrid()
{
    for (size_t i = 0; i < gridy.size(); i++)
    {
        for (size_t j = 0; j < gridx.size(); j++)
        {
            nodes.push_back(Node(gridx[j], gridy[i]));
        }
    }

    for (size_t i = 0; i < numberOfFiniteElementsy; i++)
    {
        for (size_t j = 0; j < numberOfFiniteElementsx; j++)
        {
            fes.push_back(FiniteElement(nodes[i * gridx.size() + j], nodes[i * gridx.size() + j + 1],
                nodes[(i + 1) * gridx.size() + j], nodes[(i + 1) * gridx.size() + j + 1]));
        }
    }

    globalMatrixRank = nodes.size();
    mainMatrix.resize(globalMatrixRank * globalMatrixRank);
    GMatrix.resize(globalMatrixRank * globalMatrixRank);
    MMatrix.resize(globalMatrixRank * globalMatrixRank);
    bVector.resize(globalMatrixRank);
    qVector.resize(globalMatrixRank);
}

void ForwardProblem::clear()
{
    areas.clear();
    fes.clear();
    nodes.clear();
    gridx.clear();
    gridy.clear();
    mainMatrix.clear();
    GMatrix.clear();
    MMatrix.clear();
    bVector.clear();
    qVector.clear();
    numberOfFiniteElementsx = 0;
    numberOfFiniteElementsy = 0;
    globalMatrixRank = 0;
}

std::vector<double> ForwardProblem::getMainMatrix()
{
    return mainMatrix;
}

std::vector<double> ForwardProblem::getGMatrix()
{
    return GMatrix;
}

std::vector<double> ForwardProblem::getMMatrix()
{
    return MMatrix;
}

std::vector<double> ForwardProblem::getBVector()
{
    return bVector;
}

std::vector<double> ForwardProblem::getQVector()
{
    return qVector;
}

std::vector<double> ForwardProblem::getTrueSolution()
{
    std::vector<double> returnVector(globalMatrixRank);

    for (size_t i = 0; i < nodes.size(); i++)
    {
        returnVector[i] = u(nodes[i].x, nodes[i].y);
    }

    return returnVector;
}

double ForwardProblem::getSolutionInPoint(double x, double y)
{
    int size = 4;
    for (size_t i = 0; i < fes.size(); i++)
    {
        if (x >= fes[i].leftBottom().x - eps && x <= fes[i].rightBottom().x + eps &&
            y >= fes[i].leftBottom().y - eps && y <= fes[i].leftTop().y + eps)
        {
            double sol = 0;
            std::vector<int> mapping = getGlobalPositionVector(i);

            for (size_t j = 0; j < size; j++)
            {
                sol += bifunc2d(j, Point(x, y), i) * qVector[mapping[j]];
            }
            return sol;
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double ForwardProblem::getDerivativeInPoint(double x, double y)
{
    int size = 4;
    for (size_t i = 0; i < fes.size(); i++)
    {
        if (x >= fes[i].leftBottom().x - eps && x <= fes[i].rightBottom().x + eps &&
            y >= fes[i].leftBottom().y - eps && y <= fes[i].leftTop().y + eps)
        {
            double hx = fes[i].rightBottom().x - fes[i].leftBottom().x;
            double hy = fes[i].leftTop().y - fes[i].leftBottom().y;

            double dim = hx * hy;

            double b2 = 0;
            std::vector<double> local = getLocalGMatrix(i);
            std::vector<int> mapping = getGlobalPositionVector(i);

            for (size_t j = 0; j < size; j++)
            {
                for (size_t k = 0; k < size; k++)
                {
                    b2 += local[j * size + k] * qVector[mapping[j]] * qVector[mapping[k]];
                }
            }
            b2 /= dim;
            return sqrt(b2);
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

Point ForwardProblem::getGrad(double x, double y)
{
    int size = 4;
    for (size_t i = 0; i < fes.size(); i++)
    {
        if (x >= fes[i].leftBottom().x - eps && x <= fes[i].rightBottom().x + eps &&
            y >= fes[i].leftBottom().y - eps && y <= fes[i].leftTop().y + eps)
        {
            Point p;
            std::vector<Point> fePoints = { fes[i].leftBottom(), fes[i].rightBottom(), fes[i].leftTop(), fes[i].rightTop()};
            std::vector<int> mapping = getGlobalPositionVector(i);
            for (size_t j = 0; j < size; j++)
            {
                Point grad = gradbifunc2d(j, fePoints[j], i);
                p.x += qVector[mapping[j]] * grad.x;
                p.y += qVector[mapping[j]] * grad.y;
            }
            return p;
        }
    }
    return Point();
}

void ForwardProblem::applyFirstBorder()
{
    for (size_t i = 0; i < nodes.size(); i++)
    {
        if (fabs(nodes[i].x - gridx.front()) < eps ||
            fabs(nodes[i].x - gridx.back()) < eps ||
            fabs(nodes[i].y - gridy.front()) < eps ||
            fabs(nodes[i].y - gridy.back()) < eps)
        {
            for (size_t j = 0; j < globalMatrixRank; j++)
            {
                mainMatrix[i * globalMatrixRank + j] = 0;
            }
            mainMatrix[i * globalMatrixRank + i] = 1;
            bVector[i] = u(nodes[i].x, nodes[i].y);
        }
    }
}

void ForwardProblem::applySources()
{
    int size = 4;
    for (size_t i = 0; i < fes.size(); i++)
    {
        if (source0.x >= fes[i].leftBottom().x - eps && source0.x <= fes[i].rightBottom().x + eps &&
            source0.y >= fes[i].leftBottom().y - eps && source0.y <= fes[i].leftTop().y + eps)
        {
            std::vector<int> mapping = getGlobalPositionVector(i);
            for (size_t j = 0; j < size; j++)
            {
                bVector[mapping[j]] += power * bifunc2d(j, source0, i);
            }
        }
        if (source1.x >= fes[i].leftBottom().x - eps && source1.x <= fes[i].rightBottom().x + eps &&
            source1.y >= fes[i].leftBottom().y - eps && source1.y <= fes[i].leftTop().y + eps)
        {
            std::vector<int> mapping = getGlobalPositionVector(i);
            for (size_t j = 0; j < size; j++)
            {
                bVector[mapping[j]] += -power * bifunc2d(j, source1, i);
            }
        }
    }

}

void ForwardProblem::solveProblem()
{
    calculateGMatrix();
    calculateMMatrix();
    std::transform(GMatrix.begin(), GMatrix.end(), MMatrix.begin(), mainMatrix.begin(), std::plus<double>());
    calculateBVector();
    applySources();
    applyFirstBorder();
    calculateQVector();
}

double ForwardProblem::bifunc2d(int number, Point p, int currentFE)
{
    switch (number)
    {
    case 0:
        return bifuncx(0, p.x, currentFE) * bifuncy(0, p.y, currentFE);
    case 1:
        return bifuncx(1, p.x, currentFE) * bifuncy(0, p.y, currentFE);
    case 2:
        return bifuncx(0, p.x, currentFE) * bifuncy(1, p.y, currentFE);
    case 3:
        return bifuncx(1, p.x, currentFE) * bifuncy(1, p.y, currentFE);
    default:
        break;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double ForwardProblem::bifuncx(int number, double x, int currentFE)
{
    double x0 = fes[currentFE].leftBottom().x;
    double x1 = fes[currentFE].rightBottom().x;
    double h = x1 - x0;
    if (number == 0)
    {
        return (x1 - x) / h;
    }
    else if(number == 1)
    {
        return (x - x0) / h;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double ForwardProblem::bifuncy(int number, double y, int currentFE)
{
    double y0 = fes[currentFE].leftBottom().y;
    double y1 = fes[currentFE].leftTop().y;
    double h = y1 - y0;
    if (number == 0)
    {
        return (y1 - y) / h;
    }
    else if (number == 1)
    {
        return (y - y0) / h;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

Point ForwardProblem::gradbifunc2d(int number, Point p, int currentFE)
{
    switch (number)
    {
    case 0:
        return Point(dbifuncx(0, p.x, currentFE) * bifuncy(0, p.y, currentFE), 
            bifuncx(0, p.x, currentFE) * dbifuncy(0, p.y, currentFE));
    case 1:
        return Point(dbifuncx(1, p.x, currentFE) * bifuncy(0, p.y, currentFE),
            bifuncx(1, p.x, currentFE) * dbifuncy(0, p.y, currentFE));
    case 2:
        return Point(dbifuncx(0, p.x, currentFE) * bifuncy(1, p.y, currentFE),
            bifuncx(0, p.x, currentFE) * dbifuncy(1, p.y, currentFE));
    case 3:
        return Point(dbifuncx(1, p.x, currentFE) * bifuncy(1, p.y, currentFE),
            bifuncx(1, p.x, currentFE) * dbifuncy(1, p.y, currentFE));
    default:
        break;
    }
}

double ForwardProblem::dbifuncx(int number, double x, int currentFE)
{
    double x0 = fes[currentFE].leftBottom().x;
    double x1 = fes[currentFE].rightBottom().x;
    double h = x1 - x0;
    if (number == 0)
    {
        return -1 / h;
    }
    else if (number == 1)
    {
        return 1 / h;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double ForwardProblem::dbifuncy(int number, double y, int currentFE)
{
    double y0 = fes[currentFE].leftBottom().y;
    double y1 = fes[currentFE].leftTop().y;
    double h = y1 - y0;
    if (number == 0)
    {
        return -1 / h;
    }
    else if (number == 1)
    {
        return 1 / h;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

void ForwardProblem::calculateGMatrix()
{
    for (size_t i = 0; i < fes.size(); i++)
    {
        std::vector<int> localMatrixMapping = getGlobalPositionMatrix(i);
        std::vector<double> localMatrix = getLocalGMatrix(i);
        for (size_t j = 0; j < localMatrix.size(); j++)
        {
            GMatrix[localMatrixMapping[j]] += localMatrix[j];
        }
    }
}

void ForwardProblem::calculateMMatrix()
{
    for (size_t i = 0; i < fes.size(); i++)
    {
        std::vector<int> localMatrixMapping = getGlobalPositionMatrix(i);
        std::vector<double> localMatrix = getLocalMMatrix(i);
        for (size_t j = 0; j < localMatrix.size(); j++)
        {
            MMatrix[localMatrixMapping[j]] += localMatrix[j];
        }
    }
}

void ForwardProblem::calculateBVector()
{
    for (size_t i = 0; i < fes.size(); i++)
    {
        std::vector<int> localVectorMapping = getGlobalPositionVector(i);
        std::vector<double> localVector = getLocalBVector(i);
        for (size_t j = 0; j < localVector.size(); j++)
        {
            bVector[localVectorMapping[j]] += localVector[j];
        }
    }
}

void ForwardProblem::calculateQVector()
{
    qVector = Solver::gauss(mainMatrix, bVector);
}

std::vector<double> ForwardProblem::getLocalGMatrix(int currentFE)
{
    double c1 = (getLambda(currentFE) / 6) * ((fes[currentFE].leftTop().y - fes[currentFE].leftBottom().y) /
        (fes[currentFE].rightBottom().x - fes[currentFE].leftBottom().x));
    double c2 = (getLambda(currentFE) / 6) * ((fes[currentFE].rightBottom().x - fes[currentFE].leftBottom().x) /
        (fes[currentFE].leftTop().y - fes[currentFE].leftBottom().y));
    std::vector<double> returnMatrix = { 2, -2, 1, -1, -2, 2, -1, 1, 1, -1, 2, -2, -1, 1, -2, 2 };
    std::vector<double> helpingMatrix = { 2, 1, -2, -1, 1, 2, -1, -2, -2, -1, 2, 1, -1, -2, 1, 2 };


    std::transform(returnMatrix.begin(), returnMatrix.end(), returnMatrix.begin(), [c1](double el) { return c1 * el; });
    std::transform(helpingMatrix.begin(), helpingMatrix.end(), helpingMatrix.begin(), [c2](double el) { return c2 * el; });
    std::transform(returnMatrix.begin(), returnMatrix.end(), helpingMatrix.begin(), returnMatrix.begin(), std::plus<double>());

    return returnMatrix;
}

std::vector<double> ForwardProblem::getLocalCMatrix(int currentFE)
{
    double c1 = ((fes[currentFE].rightBottom().x - fes[currentFE].leftBottom().x) *
        (fes[currentFE].leftTop().y - fes[currentFE].leftBottom().y)) / 36;
    std::vector<double> returnMatrix = { 4, 2, 2, 1, 2, 4, 1, 2, 2, 1, 4, 2, 1, 2, 2, 4 };
    std::transform(returnMatrix.begin(), returnMatrix.end(), returnMatrix.begin(), [c1](double el) { return c1 * el; });
    return returnMatrix;
}

std::vector<double> ForwardProblem::getLocalMMatrix(int currentFE)
{
    std::vector<double> returnMatrix = getLocalCMatrix(currentFE);
    std::transform(returnMatrix.begin(), returnMatrix.end(), returnMatrix.begin(), [&](double el) { return getGamma(currentFE) * el; });
    return returnMatrix;
}

std::vector<double> ForwardProblem::getLocalBVector(int currentFE)
{
    int size = 4;
    std::vector<double> returnVector(size);
    std::vector<double> cMatrix = getLocalCMatrix(currentFE);
    std::vector<double> fVector = { f(fes[currentFE].leftBottom().x, fes[currentFE].leftBottom().y, getLambda(currentFE), getGamma(currentFE)),
    f(fes[currentFE].rightBottom().x, fes[currentFE].rightBottom().y, getLambda(currentFE), getGamma(currentFE)),
    f(fes[currentFE].leftTop().x, fes[currentFE].leftTop().y, getLambda(currentFE), getGamma(currentFE)),
    f(fes[currentFE].rightTop().x, fes[currentFE].rightTop().y, getLambda(currentFE), getGamma(currentFE)) };

    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            returnVector[i] += cMatrix[i * size + j] * fVector[j];
        }
    }

    return returnVector;
}

std::vector<int> ForwardProblem::getGlobalPositionMatrix(int currentFE)
{
    int size = 4;
    std::vector<int> mappingNodes;
    int firstNode = (currentFE + currentFE / numberOfFiniteElementsx);
    int secondNode = (currentFE + currentFE / numberOfFiniteElementsx + 1);
    int thirdNode = (currentFE + numberOfFiniteElementsx + currentFE / numberOfFiniteElementsx + 1);
    int forthNode = (currentFE + numberOfFiniteElementsx + currentFE / numberOfFiniteElementsx + 1 + 1);

    mappingNodes = { firstNode, secondNode, thirdNode, forthNode };

    std::vector<int> returnIndices(size * size);

    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            returnIndices[i * size + j] = mappingNodes[i] * globalMatrixRank + mappingNodes[j];
        }
    }

    return returnIndices;
}

std::vector<int> ForwardProblem::getGlobalPositionVector(int currentFE)
{
    int size = 4;
    std::vector<int> mappingNodes;
    int firstNode = (currentFE + currentFE / numberOfFiniteElementsx);
    int secondNode = (currentFE + currentFE / numberOfFiniteElementsx + 1);
    int thirdNode = (currentFE + numberOfFiniteElementsx + currentFE / numberOfFiniteElementsx + 1);
    int forthNode = (currentFE + numberOfFiniteElementsx + currentFE / numberOfFiniteElementsx + 1 + 1);

    mappingNodes = { firstNode, secondNode, thirdNode, forthNode };

    return mappingNodes;
}

double ForwardProblem::getLambda(int currentFE)
{
    for (size_t i = 0; i < areas.size(); i++)
    {
        if (areas[i].contains(fes[currentFE]))
        {
            return areas[i].getLambda();
        }
    }
    return 0.0;
}

double ForwardProblem::getGamma(int currentFE)
{
    for (size_t i = 0; i < areas.size(); i++)
    {
        if (areas[i].contains(fes[currentFE]))
        {
            return areas[i].getGamma();
        }
    }
    return 0.0;
}

double ForwardProblem::getSigma(int currentFE)
{
    for (size_t i = 0; i < areas.size(); i++)
    {
        if (areas[i].contains(fes[currentFE]))
        {
            return areas[i].getSigma();
        }
    }
    return 0.0;
}
