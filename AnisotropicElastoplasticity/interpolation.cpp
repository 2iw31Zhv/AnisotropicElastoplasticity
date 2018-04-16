#include "interpolation.h"
#include "interpolation.h"
#include <cmath>
#include <functional>

using namespace std;
using namespace Eigen;

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


double hue2rgb(double p, double q, double t)
{
	if (t < 0) t += 1.0;
	if (t > 1.0) t -= 1.0;
	if (t < 1.0 / 6.0) return p + (q - p) * 6 * t;
	if (t < 0.5) return q;
	if (t < 2.0 / 3.0) return p + (q - p)*(2.0 / 3.0 - t) * 6;
	return p;
}

Vector3d hslToRgb(double h, double s, double l)
{
	Vector3d result(l, l, l);
	if (s == 0)
		return result;
	double q = (l < 0.5 ? l*(1.0 + s) : l + s - l*s);
	double p = 2.0*l - q;
	result[0] = hue2rgb(p, q, h + 1.0 / 3.0);
	result[1] = hue2rgb(p, q, h);
	result[2] = hue2rgb(p, q, h - 1.0 / 3.0);
	return result;
}