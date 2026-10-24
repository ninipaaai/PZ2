#include <iostream>
#include <cmath>
#include <stdexcept>
#include "Cubic_Interpolation_Spline_1D.h"

namespace Com_Methods
{
	void Cubic_Interpolation_Spline_1D::Update_Spline(const std::vector<Point> &Points, 
		                                              const std::vector<double> &F_Value)
	{
		this->Points.clear();
		for (auto & x : Points) this->Points.push_back(x);
		int n = Points.size() - 1;
		double h1, h2;
		a.resize(n); b.resize(n); c.resize(n); d.resize(n);
		std::vector<double> f(n - 1);

		for (int i = 0; i < n - 1; i++)
		{
			h1 = Points[i + 1].x() - Points[i].x();
			h2 = Points[i + 2].x() - Points[i + 1].x();
			b[i] = 2 * (h1 + h2);
			a[i + 1] = h1;
			d[i] = h2;
			f[i] = 3.0 * ((F_Value[i + 2] - F_Value[i + 1]) / h2 - (F_Value[i + 1] - F_Value[i]) / h1);
		}
		
		for (int j = 1; j < n - 1; j++)
		{
			b[j] -= a[j] / b[j - 1] * d[j - 1];
			f[j] -= a[j] / b[j - 1] * f[j - 1]; 
		}

		c[n - 1] = f[n - 2] / b[n - 2];
		for (int j = n - 2; j > 0; j--)
		{
			c[j] = (f[j - 1] - c[j + 1] * d[j - 1]) / b[j - 1];
		}

		c[0] = 0.0;

		for (int i = 0; i < n - 1; i++)
		{
			h1 = Points[i + 1].x() - Points[i].x();
			a[i] = F_Value[i];
			b[i] = (F_Value[i + 1] - F_Value[i]) / h1 - (c[i + 1] + 2.0 * c[i]) * h1 / 3.0;
			d[i] = (c[i + 1] - c[i]) / h1 / 3.0;
		}

		h1 = Points[n].x() - Points[n - 1].x();
		a[n - 1] = F_Value[n - 1];
		b[n - 1] = (F_Value[n] - F_Value[n - 1]) / h1 - 2.0 * c[n - 1] * h1 / 3.0;
		d[n - 1] = -c[n - 1] / h1 / 3.0;
	}

	
	void Cubic_Interpolation_Spline_1D::Get_Value(const Point &P, double * data)const
	{
		double eps = 1e-7;
		int n = Points.size() - 1;
		for (int i = 0; i < n; i++)
		{
			if (P.x() > Points[i].x() && P.x() < Points[i + 1].x() || fabs(P.x() - Points[i].x()) < eps || fabs(P.x() - Points[i + 1].x()) < eps)
			{
				double diff = (P.x() - Points[i].x());
				data[0] = a[i] + b[i] * diff + c[i] * pow(diff, 2) + d[i] * pow(diff, 3); 
				data[1] = b[i] + 2.0 * c[i] * diff + 3.0 * d[i] * pow(diff, 2);
				data[2] = 2.0 * c[i] + 6.0 * d[i] * diff;
				return;
			}
		}
	}
}

