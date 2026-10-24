#pragma once

#include "Point.h"
#include <vector> 

//Формирование одномерной сетки
namespace Com_Methods
{
    std::vector<Point> Grid_1D(double a, double b, int count, double r = 1.0);
}