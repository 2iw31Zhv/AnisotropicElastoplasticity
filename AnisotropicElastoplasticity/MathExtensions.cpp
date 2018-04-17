#include <cmath>
#include <iostream>
#include "MathExtensions.h"


namespace MathExtensions
{
	double acos(double x)
	{
		if (x > 1.0)
		{
			std::cerr << "acos: x > 1.0!\n";
			exit(1);
		}
		else if (x < -1.0)
		{
			std::cerr << "acos: x < -1.0!\n";
			exit(1);
		}

		return ::acos(x);
	}
	double dacos(double x)
	{
		if (fabs(1.0 - x * x) < 1e-5)
		{
			std::cerr << "dacos fabs(1 - x * x) < 1e-5!\n";
			exit(1);
		}
		return -1.0 / sqrt(1.0 - x * x);
	}
	double ddacos(double x)
	{
		if (fabs(1.0 - x * x) < 1e-5)
		{
			std::cerr << "ddacos fabs(1 - x * x) < 1e-5!\n";
			exit(1);
		}
		return - x / pow(1.0 - x * x, 1.5);
	}

	double asin(double x)
	{
		if (x < -1.0)
		{
			std::cerr << "asin: x < -1.0!\n";
		}
		else if (x > 1.0)
		{
			std::cerr << "asin: x > 1.0!\n";
		}

		return ::asin(x);
	}
	double dasin(double x)
	{
		if (fabs(1.0 - x * x) < 1e-5 || x * x > 1.0)
		{
			std::cerr << "dasin fabs(1 - x * x) < 1e-5!\n";
			exit(1);
		}
		return 1.0 / sqrt(1 - x * x);
	}
	double ddasin(double x)
	{
		if (fabs(1.0 - x * x) < 1e-5)
		{
			std::cerr << "ddasin fabs(1 - x * x) < 1e-5!\n";
			exit(1);
		}
		return x / pow(1.0 - x * x, 1.5);
	}
}