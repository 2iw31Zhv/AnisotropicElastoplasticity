#include "interpolation.h"
#include "interpolation.h"
#include <cmath>
#include <functional>

using namespace std;

double cubic_B_spline(double x)
{
	double abs_x = abs(x);
	return (abs_x >= 2.0) ? 0.0 : 
		(abs_x >= 1.0 ? -1.0 / 6.0 * abs_x * abs_x * abs_x
			+ abs_x * abs_x - 2.0 * abs_x + 4.0 / 3.0 : 
			0.5 * abs_x * abs_x * abs_x - abs_x * abs_x + 2.0 / 3.0);
}

double Dcubic_B_spline(double x)
{
	return (x >= 2.0) ? 0.0 :
		(
			x >= 1.0 ? -0.5 * x * x + 2.0 * x - 2.0 :
			(
				x >= 0.0 ? 1.5 * x * x - 2.0 * x :
				(
					x >= -1.0 ? -1.5 * x * x - 2.0 * x :
					(
						x >= -2.0 ? 0.5 * x * x + 2.0 * x + 2.0 : 0.0
						)
					)
				)
			);
}

double clamp(double x, double lowerBound, double higherBound)
{
	if (x > higherBound)
	{
		return higherBound;
	}
	else if (x < lowerBound)
	{
		return lowerBound;
	}
	else
	{
		return x;
	}
}
