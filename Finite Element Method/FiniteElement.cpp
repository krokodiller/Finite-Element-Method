#include "FiniteElement.h"
constexpr double eps = 0.00001;
Point::Point() : x(0), y(0)
{

}

Point::Point(double x, double y) : x(x), y(y)
{

}

bool Point::operator==(const Point& other)
{
    return this->x == other.x && this->y == other.y;
}

FiniteElement::FiniteElement() : nodes(4)
{

}

FiniteElement::FiniteElement(Node leftbottom, Node rightbottom, Node lefttop, Node righttop) : nodes(4)
{
    nodes[0] = leftbottom;
    nodes[1] = rightbottom;
    nodes[2] = lefttop;
    nodes[3] = righttop;
}

Node FiniteElement::leftBottom()
{
    return nodes[Points::leftbottom];
}

Node FiniteElement::rightBottom()
{
    return nodes[Points::rightbottom];
}

Node FiniteElement::leftTop()
{
    return nodes[Points::lefttop];
}

Node FiniteElement::rightTop()
{
    return nodes[Points::righttop];
}

TriangleFiniteElement::TriangleFiniteElement() : nodes(3)
{

}

TriangleFiniteElement::TriangleFiniteElement(Node first, Node second, Node third) : nodes(3)
{
    nodes[0] = first;
    nodes[1] = second;
    nodes[2] = third;
}

Node TriangleFiniteElement::first()
{
    return nodes[0];
}

Node TriangleFiniteElement::second()
{
    return nodes[1];
}

Node TriangleFiniteElement::third()
{
    return nodes[2];
}

Area::Area() : points(4), lambda(0), gamma(0)
{
    
}

Area::Area(double xleftbottom, double yleftbottom, double xrighttop, double yrighttop) :
    points(4), lambda(0), gamma(0)
{
    points[0] = Point(xleftbottom, yleftbottom);
    points[1] = Point(xrighttop, yleftbottom);
    points[2] = Point(xleftbottom, yrighttop);
    points[3] = Point(xrighttop, yrighttop);
}

Area::Area(double xleftbottom, double yleftbottom, double xrighttop, double yrighttop, double lambda, double gamma) :
     points(4), lambda(lambda), gamma(gamma)
{
    points[0] = Point(xleftbottom, yleftbottom);
    points[1] = Point(xrighttop, yleftbottom);
    points[2] = Point(xleftbottom, yrighttop);
    points[3] = Point(xrighttop, yrighttop);
}

Area::Area(double xleftbottom, double yleftbottom, double xrighttop, double yrighttop, double sigma) :
    points(4), sigma(sigma), gamma(0)
{
    points[0] = Point(xleftbottom, yleftbottom);
    points[1] = Point(xrighttop, yleftbottom);
    points[2] = Point(xleftbottom, yrighttop);
    points[3] = Point(xrighttop, yrighttop);
}

bool Area::contains(double x, double y)
{
    if (x >= points[leftbottom].x - eps  && x <= points[rightbottom].x + eps && 
        y >= points[leftbottom].y - eps && y <= points[lefttop].y + eps)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Area::contains(Point p)
{
    if (p.x >= points[leftbottom].x - eps && p.x <= points[rightbottom].x + eps &&
        p.y >= points[leftbottom].y - eps && p.y <= points[lefttop].y + eps)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Area::contains(FiniteElement fe)
{
    if (fe.leftBottom().x >= points[leftbottom].x - eps && fe.rightBottom().x <= points[rightbottom].x + eps &&
        fe.leftBottom().y >= points[leftbottom].y - eps && fe.leftTop().y <= points[lefttop].y + eps)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Area::setLambda(double lambda)
{
    this->lambda = lambda;
}

void Area::setGamma(double gamma)
{
    this->gamma = gamma;
}

void Area::setSigma(double sigma)
{
    this->sigma = sigma;
}

double Area::getLambda()
{
    return lambda;
}

double Area::getGamma()
{
    return gamma;
}

double Area::getSigma()
{
    return sigma;
}
