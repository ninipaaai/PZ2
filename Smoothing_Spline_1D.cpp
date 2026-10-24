#include <functional>
#include "Smoothing_Spline_1D.h"
#include <iostream>
#include <cmath>

namespace Com_Methods
{
	Smoothing_Spline_1D::Smoothing_Spline_1D(const double &p)
	{
		this->p = p;
	}

	void Smoothing_Spline_1D::Transition_To_Master_Element(int segment, const double &x, double &ksi) const
	{
		ksi = 2.0 * (x - Points[segment].x()) / (Points[segment + 1].x() - Points[segment].x()) - 1.0;
	}

	double Smoothing_Spline_1D::Basis_Function(int Number, const double &Ksi) const
	{
		switch (Number) {
			case 1:  {return 0.5*(1 - Ksi); break; }
			case 2:  {return 0.5*(1 + Ksi); break; }
			default: {throw std::runtime_error("Ошибка базисной функции"); break; }
		}
	}

	double Smoothing_Spline_1D::Der_Basis_Function(int number, const double &ksi) const
	{
		switch (number) {
			case 1:  {return -0.5; break; }
			case 2:  {return  0.5; break; }
			default: {throw std::runtime_error("Ошибка базисной функции"); break; }
		}
	}

	void Smoothing_Spline_1D::Update_Spline(const std::vector<Point>  & Points,
											const std::vector<double> & F_Value)
	{
		this->Points.clear();
		for (auto & elem : Points) this->Points.push_back(elem);
		int segments = Points.size() - 1;
		alpha.resize(segments + 1);

		std::vector <double> a, b, c;
		a.resize(segments + 1); b.resize(segments + 1); c.resize(segments + 1);

		std::function<void(int segments, const Point &point, const double &F_Val, const double &w)> 
		Assembling = [&](int i, const Point &point, const double &F_Val, const double &w) {
			double X = point.x(), Ksi;
		    Transition_To_Master_Element(i, X, Ksi);
			double f1 = Basis_Function(1, Ksi);
			double f2 = Basis_Function(2, Ksi);

			b[i]     += (1.0 - p) * w * f1 * f1;
			b[i + 1] += (1.0 - p) * w * f2 * f2;
			a[i + 1] += (1.0 - p) * w * f1 * f2;
			c[i]	 += (1.0 - p) * w * f2 * f1;
			alpha[i]     += (1.0 - p) * w * f1 * F_Val;
			alpha[i + 1] += (1.0 - p) * w * f2 * F_Val;
		};

		for (int i = 0; i < segments; i++) {
			double W = 1;
			Assembling(i, this->Points[i], F_Value[i], W);
			Assembling(i, this->Points[i + 1], F_Value[i + 1], W);

			double h = Points[i + 1].x() - Points[i].x();
			b[i]	 += 1.0 / h * p;
			b[i + 1] += 1.0 / h * p;
			a[i + 1] -= 1.0 / h * p;
			c[i]	 -= 1.0 / h * p;
		}

		for (int j = 1; j < segments + 1; j++) {
			b[j] -= a[j] / b[j - 1] * c[j - 1];         
			alpha[j] -= a[j] / b[j - 1] * alpha[j - 1]; 
		}

		alpha[segments] /= b[segments];
		for (int j = segments - 1; j >= 0; j--) alpha[j] = (alpha[j] - alpha[j + 1] * c[j]) / b[j];
	}

	void Smoothing_Spline_1D::Get_Value(const Point &point, double * Res)const
	{
		double eps = 1e-7;
		int segments = Points.size() - 1;
		double X = point.x();

		for (int i = 0; i < segments; i++) {
			if (X > Points[i].x() && X < Points[i + 1].x() ||
				fabs(X - Points[i].x()) < eps ||
				fabs(X - Points[i + 1].x()) < eps)
			{
				double h = Points[i + 1].x() - Points[i].x();
				double Ksi;
				Transition_To_Master_Element(i, X, Ksi);
				Res[0] = alpha[i]* Basis_Function(1, Ksi) + alpha[i + 1] * Basis_Function(2, Ksi);
				Res[1] = (alpha[i]* Der_Basis_Function(1, Ksi) + alpha[i + 1] * Der_Basis_Function(2, Ksi)) * 2.0 / h;
				Res[2] = 0.0;
				return;
			}
		}
		throw std::runtime_error("The point is not found in the segments...");
	}
}
