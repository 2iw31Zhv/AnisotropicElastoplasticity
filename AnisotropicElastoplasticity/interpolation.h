#pragma once

#include <Eigen/Core>

double cubic_B_spline(double x);

double Dcubic_B_spline(double x);

double clamp(double x, double lowerBound, double higherBound);

double hue2rgb(double p, double q, double t);

Eigen::Vector3d hslToRgb(double h, double s, double l);

namespace Spectrum
{
	const Eigen::Vector3d red(1.0, 0.0, 0.0);
	const Eigen::Vector3d orange(1.0, 0.647, 0.0);
	const Eigen::Vector3d yellow(1.0, 1.0, 0.0);
	const Eigen::Vector3d green(0.0, 1.0, 0.0);
	const Eigen::Vector3d blue(0.0, 0.498, 1.0);
	const Eigen::Vector3d indigo(0.0, 0.0, 1.0);
	const Eigen::Vector3d purple(0.545, 0.0, 1.0);
}