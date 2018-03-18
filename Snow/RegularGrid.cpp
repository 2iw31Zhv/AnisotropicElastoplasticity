#include "RegularGrid.h"
#include <iostream>
#include <igl/viewer/Viewer.h>

using namespace std;
using namespace Eigen;
using namespace igl;

static double hue2rgb(double p, double q, double t)
{
	if (t < 0) t += 1.0;
	if (t > 1.0) t -= 1.0;
	if (t < 1.0 / 6.0) return p + (q - p) * 6 * t;
	if (t < 0.5) return q;
	if (t < 2.0 / 3.0) return p + (q - p)*(2.0 / 3.0 - t) * 6;
	return p;
}

static Vector3d hslToRgb(double h, double s, double l)
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

void RegularGrid::initializeRenderingData_()
{
	colors_.resize(gridNumber() * 12, 3);
	colors_.col(2).setOnes();
	points_1_.resize(gridNumber() * 12, 3);
	points_2_.resize(gridNumber() * 12, 3);

	auto VTX = [&](int i, int j, int k)
	{
		return Vector3d(i * h_[0], j * h_[1], k * h_[2]);
	};

	for (int k = 0; k < resolution_[2]; ++k)
	{
		for (int j = 0; j < resolution_[1]; ++j)
		{
			for (int i = 0; i < resolution_[0]; ++i)
			{

#define VTX(s0, s1, s2) \
Vector3d v##s0##s1##s2(minBound_[0] + (i + s0 - 0.5) * h_[0], minBound_[1] + (j + s1 - 0.5) * h_[1], minBound_[2] + (k + s2 - 0.5) * h_[2]);

				VTX(0, 0, 0);
				VTX(0, 0, 1);
				VTX(0, 1, 0);
				VTX(0, 1, 1);
				VTX(1, 0, 0);
				VTX(1, 0, 1);
				VTX(1, 1, 0);
				VTX(1, 1, 1);
#undef VTX

				int index = toIndex(i, j, k);

				points_1_.row(12 * index + 0) = v000;
				points_1_.row(12 * index + 1) = v100;
				points_1_.row(12 * index + 2) = v110;
				points_1_.row(12 * index + 3) = v010;

				points_1_.row(12 * index + 4) = v000;
				points_1_.row(12 * index + 5) = v100;
				points_1_.row(12 * index + 6) = v110;
				points_1_.row(12 * index + 7) = v010;

				points_1_.row(12 * index + 8) = v001;
				points_1_.row(12 * index + 9) = v101;
				points_1_.row(12 * index + 10) = v111;
				points_1_.row(12 * index + 11) = v011;

				points_2_.row(12 * index + 0) = v100;
				points_2_.row(12 * index + 1) = v110;
				points_2_.row(12 * index + 2) = v010;
				points_2_.row(12 * index + 3) = v000;

				points_2_.row(12 * index + 4) = v001;
				points_2_.row(12 * index + 5) = v101;
				points_2_.row(12 * index + 6) = v111;
				points_2_.row(12 * index + 7) = v011;

				points_2_.row(12 * index + 8) = v101;
				points_2_.row(12 * index + 9) = v111;
				points_2_.row(12 * index + 10) = v011;
				points_2_.row(12 * index + 11) = v001;

			}
		}
	}
}

void RegularGrid::recomputeColors()
{
	double massMax = 1e-2;
	for (int k = 0; k < resolution_[2]; ++k)
	{
		for (int j = 0; j < resolution_[1]; ++j)
		{
			for (int i = 0; i < resolution_[0]; ++i)
			{

				int index = toIndex(i, j, k);

				double colorVal = min(masses[index] / massMax, 1.0);
				double h = (1.0 - colorVal) * 240.0 / 360.0;
				double s = 100.0 / 100.0;
				double l = 50.0 / 100.0;

				Vector3d color = hslToRgb(h, s, l);

				colors_.row(12 * index + 0) = color;
				colors_.row(12 * index + 1) = color;
				colors_.row(12 * index + 2) = color;
				colors_.row(12 * index + 3) = color;
				colors_.row(12 * index + 4) = color;
				colors_.row(12 * index + 5) = color;
				colors_.row(12 * index + 6) = color;
				colors_.row(12 * index + 7) = color;
				colors_.row(12 * index + 8) = color;
				colors_.row(12 * index + 9) = color;
				colors_.row(12 * index + 10) = color;
				colors_.row(12 * index + 11) = color;
			}
		}
	}
}

RegularGrid::RegularGrid(
	const Eigen::VectorXd & minBound, 
	const Eigen::VectorXd & maxBound, 
	const Eigen::Vector3i & resolution):
	minBound_(minBound),
	maxBound_(maxBound),
	resolution_(resolution)
{
	if (maxBound[0] < minBound[0] ||
		maxBound[1] < minBound[1] ||
		maxBound[2] < maxBound[2])
	{
		std::cerr << "[ERROR] bad parameter: "
			<< __FILE__ << ", "
			<< __FUNCTION__ << ", "
			<< __LINE__ << ": "
			<< "maxBound must be bigger than minBound!\n";
		while (1);
	}

	h_[0] = (maxBound_[0] - minBound_[0]) / resolution_[0];
	h_[1] = (maxBound_[1] - minBound_[1]) / resolution_[1];
	h_[2] = (maxBound_[2] - minBound_[2]) / resolution_[2];

	masses.resize(gridNumber());
	forces.resize(gridNumber(), 3);
	velocities.resize(gridNumber(), 3);

	initializeRenderingData_();
}

int RegularGrid::toIndex(int i, int j, int k) const
{
	return k * resolution_[0] * resolution_[1]
		+ j * resolution_[0] + i;
}

std::tuple<int, int, int> RegularGrid::toCoordinate(int index) const
{
	int k = index / (resolution_[0] * resolution_[1]);
	int j = (index % (resolution_[0] * resolution_[1])) / resolution_[0];
	int i = index - k * resolution_[0] * resolution_[1]
		- j * resolution_[0];
	return make_tuple(i, j, k);
}

void RegularGrid::updateViewer()
{
	using namespace viewer;
	viewer_->data.add_edges(points_1_, points_2_, colors_);
}
