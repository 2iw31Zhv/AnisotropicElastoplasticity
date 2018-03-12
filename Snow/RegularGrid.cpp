#include "RegularGrid.h"
#include <iostream>
#include <igl/viewer/Viewer.h>

using namespace std;
using namespace Eigen;
using namespace igl;

int RegularGrid::toIndex_(int i, int j, int k) const
{
	return k * resolution_[0] * resolution_[1]
		+ j * resolution_[0] + i;
}

void RegularGrid::initializeRenderingData_()
{
	colors_.resize(gridNumber(), 3);
	colors_.col(2).setOnes();
	points_1_.resize(gridNumber() * 12, 3);
	points_2_.resize(gridNumber() * 12, 3);

	auto VTX = [&](int i, int j, int k)
	{
		return Vector3d(i * h_[0], j * h_[1], k * h_[2]);
	};

	for (int i = 0; i < resolution_[0]; ++i)
	{
		for (int j = 0; j < resolution_[1]; ++j)
		{
			for (int k = 0; k < resolution_[2]; ++k)
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

				int index = toIndex_(i, j, k);

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

	data_.resize(resolution_[0] * resolution_[1] * resolution_[2]);

	initializeRenderingData_();
}

const GridValue & RegularGrid::at(int i, int j, int k) const
{
	return data_[toIndex_(i, j, k)];
}

void RegularGrid::set(int i, int j, int k, const GridValue & gridvalue)
{
	data_[toIndex_(i, j, k)] = gridvalue;
}



void RegularGrid::updateViewer()
{
	using namespace viewer;
	viewer_->data.add_edges(points_1_, points_2_, colors_);
}
