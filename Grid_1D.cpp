#include "Grid_1D.h"
#include <cmath>

namespace Com_Methods
{
    std::vector<Point> Grid_1D(double a, double b, int count, double r)
    {
        std::vector<Point> result;
        result.reserve(count);

        if (r == 1.0)
        {
            //равномерная сетка: x_i = a + i*h
            double h = (b - a) / (count - 1);
            for (int i = 0; i < count; i++)
                result.emplace_back(a + i * h, 0.0, 0.0);
        }
        else
        {
            //неравномерная сетка: шаги образуют геометрическую прогрессию
            double sum = 0.0;
            for (int i = 0; i < count - 1; i++)
                sum += std::pow(r, i);
            double h0 = (b - a) / sum;

            double x = a;
            result.emplace_back(x, 0.0, 0.0);
            for (int i = 1; i < count; i++)
            {
                x += h0 * std::pow(r, i - 1);
                result.emplace_back(x, 0.0, 0.0);
            }
        }

        return result;
    }
}