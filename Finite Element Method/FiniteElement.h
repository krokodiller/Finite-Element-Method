#pragma once
#include <vector>


class Point
{
public:
    Point();
    Point(double x, double y);
    bool operator==(const Point &other);

    double x, y;
};

using Node = Point;
class FiniteElement
{
public:
    FiniteElement();
    FiniteElement(Node leftbottom, Node rightbottom, Node lefttop, Point righttop);
    Node leftBottom();
    Node rightBottom();
    Node leftTop();
    Node rightTop();

private:
    std::vector<Node> nodes;
};

using Node = Point;
class TriangleFiniteElement
{
public:
    TriangleFiniteElement();
    TriangleFiniteElement(Node first, Node second, Node third);
    Node first();
    Node second();
    Node third();
private:
    std::vector<Node> nodes;

};

class Area
{
public:
    Area();
    Area(double xleftbottom, double yleftbottom, double xrighttop, double yrighttop);
    Area(double xleftbottom, double yleftbottom, double xrighttop, double yrighttop, double lambda, double gamma);
    Area(double xleftbottom, double yleftbottom, double xrighttop, double yrighttop, double sigma);

    bool contains(double x, double y);
    bool contains(Point p);
    bool contains(FiniteElement fe);

    void setLambda(double lambda);
    void setGamma(double gamma);
    void setSigma(double sigma);

    double getLambda();
    double getGamma();
    double getSigma();

private:
    std::vector<Point> points;
    union
    {
        double lambda;
        double sigma;
    };

    double gamma;
};


enum Points
{
    leftbottom,
    rightbottom,
    lefttop,
    righttop
};