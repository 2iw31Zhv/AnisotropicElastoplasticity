#include "LinearAlgebra.h"
#include <cmath>

using namespace Eigen;

Eigen::Matrix3d crossMatrix(const Eigen::Vector3d & e)
{
	Matrix3d skewMatrix;
	skewMatrix(0, 0) = 0;
	skewMatrix(0, 1) = -e.z();
	skewMatrix(0, 2) = e.y();
	
	skewMatrix(1, 0) = e.z();
	skewMatrix(1, 1) = 0;
	skewMatrix(1, 2) = -e.x();

	skewMatrix(2, 0) = -e.y();
	skewMatrix(2, 1) = e.x();
	skewMatrix(2, 2) = 0;

	return skewMatrix;
}

