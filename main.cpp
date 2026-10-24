#include <iostream>
#include <windows.h>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "Smoothing_Spline_1D.h"
#include "Cubic_Interpolation_Spline_1D.h"
#include "Grid_1D.h"
using namespace std;

const double EPS=1e-12;

//Моделирование стандартного нормального распределения c помощью метода Бокса-Мюллера
double boxMuller(std::mt19937& gen) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double u1 = uniform(gen); //равномерно распределены
    double u2 = uniform(gen);

    if (u1 < EPS) u1 = EPS;  //защита от log(0)

    double r = std::sqrt(-2.0 * std::log(u1));
    double theta = 2.0 * M_PI * u2;             

    return r * std::cos(theta);
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8); 
    
    //Вариант 13
    const int N = 1784; //объём выборки
    const double mean = 1.01;
    const double scale = 3.85;
    const double a = 0; const double b = 100;
    const int w = 15;
    cout << "Практическое задание №2. Обухов Назар ПМИ-41" << endl;

    vector<Com_Methods::Point> grid = Com_Methods::Grid_1D(0, 100, N); //формирование сетки
    std::random_device rd;
    std::mt19937 gen(42);
    double x;
    vector<double> F_val; F_val.reserve(N); //значения функции
    for (int i = 0; i<N; i++) {
        x = boxMuller(gen)*scale+mean;
        F_val.push_back(x);
    }

    Com_Methods::Smoothing_Spline_1D spline1(0.5); //задаём сплайны
    spline1.Update_Spline(grid, F_val);
    Com_Methods::Cubic_Interpolation_Spline_1D spline2;
    spline2.Update_Spline(grid, F_val);
    vector<double> spline1_value; vector<double> spline2_value;

    cout << "         x    Сглаживающий сплайн    Интерполяционный сплайн" << endl;
    double s1[3], s2[3];
    for (int i=0; i<N; i++) {
        spline1.Get_Value(grid[i], s1);
        spline2.Get_Value(grid[i], s2);
        spline1_value.push_back(s1[0]); spline2_value.push_back(s2[0]);
        if (i<=10) cout << setw(w) << grid[i].x() << setw(w) << spline1_value[i] << setw(w) << spline2_value[i] << endl;
    }


    return 0;
}