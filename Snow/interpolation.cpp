#include "interpolation.h"
#include "interpolation.h"
#include <cmath>

double cubic_B_spline(double x)
{
	double abs_x = abs(x);
	if (abs_x >= 2.0)
	{
		return 0.0;
	}
	else if (1.0 <= abs_x && abs_x < 2)
	{
		return -1.0 / 6.0 * abs_x * abs_x * abs_x 
			+ abs_x * abs_x - 2.0 * abs_x + 4.0 / 3.0;
	}
	else
	{
		return 0.5 * abs_x * abs_x * abs_x - abs_x * abs_x + 2.0 / 3.0;
	}
}
