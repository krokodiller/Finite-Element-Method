#include "FiniteElement.h"
#include "ForwardProblem.h"
#include <fstream>
#include <iostream>
constexpr double size = 32;

double deriv(double x, double y)
{
    double dx = 3 * x * x * y * y * y;
    double dy = x * x * x * 3 * y * y;
    double ans = dx * dx + dy * dy;
    return sqrt(ans);
}

int main()
{
    ForwardProblem fw;
    fw.setGridx(0, 4, size, 1);
    fw.setGridy(0, 4, size, 1);
    fw.addArea(0, 0, 4, 4, 1, 1);

    fw.assembleGrid();
    fw.solveProblem();


    std::ofstream s("femsolution.txt");
    std::ofstream d("femderivatives.txt");
    double h = 0.5;
    for (double i = 0; i < 4; i += 4 / size)
    {
        for (double j = 0; j < 4; j += 4 / size)
        {
            s << i << " " << j << " " << fw.getSolutionInPoint(i, j) << " " << 1 << "\n";
        }
    }
    //for (double i = 0, j = 0; i < 4; i+=0.5, j+=0.5)
    //{
    //    s << i << " " << j << " " << fw.getSolutionInPoint(i, j) << " " << 1 << "\n";
    //    d << i << " " << j << " " << fw.getDerivativeInPoint(i, j) << "\n";
    //}
    s.close();

    std::vector<double> x{ 0.1, 0.1, 1.7, 2.5, 2.5, 1.5 };
    std::vector<double> y{ 0.1, 3.9, 0.1, 2.5, 3.9, 2.5 };

    for (size_t i = 0; i < x.size(); i++)
    {
        d.precision(10);
        d << x[i] << " " << y[i] << " " << fw.getDerivativeInPoint(x[i], y[i]) << "\n";
    }
    d.close();

    std::ofstream q("solpoints.txt");
    for (size_t i = 0; i < x.size(); i++)
    {
        q.precision(10);
        q << fw.getSolutionInPoint(x[i], y[i]) << "\n";
    }
    q.close();
    std::ofstream w("derpoints.txt");
    for (size_t i = 0; i < x.size(); i++)
    {
        w.precision(10);
        w << fw.getDerivativeInPoint(x[i], y[i]) << "\n";
    }

    std::ofstream e("truederpoints.txt");
    for (size_t i = 0; i < x.size(); i++)
    {
        e.precision(10);
        e << deriv(x[i], y[i]) << "\n";
    }

    std::ofstream r("truesolpoints.txt");
    for (size_t i = 0; i < x.size(); i++)
    {
        r.precision(10);
        r << fw.u(x[i], y[i]) << "\n";
    }

    return 0;
}