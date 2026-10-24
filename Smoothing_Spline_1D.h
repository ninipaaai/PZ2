#pragma once
#ifndef Smoothing_Spline_1D_h
#define Smoothing_Spline_1D_h

#include "Spline.h"

namespace Com_Methods
{
	class Smoothing_Spline_1D : public Spline
	{
	private:
		double p; 	//параметр сглаживания
		std::vector<Point> Points; //сетка
		std::vector<double> alpha; //разложение по базису
		void Transition_To_Master_Element(int segment, const double &x, double &ksi) const; //Переход к отрезку [-1; 1]
		double Basis_Function(int number, const double &ksi) const; //локальные базисные функции
		double Der_Basis_Function(int number, const double & ksi) const; //Производные от базисных функций
	public:
		Smoothing_Spline_1D(const double &p); //конструктор
		void Update_Spline(const std::vector<Point> &Points, const std::vector<double> &F_Value) override;
		void Get_Value(const Point &P, double * Res)const override; //значение в точке
	};
}

#endif