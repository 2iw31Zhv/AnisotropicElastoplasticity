#include "geometry.h"
#include <Eigen/Dense>

using namespace std;

double segmentLength(const Eigen::Vector3d & v1, const Eigen::Vector3d & v2)
{
	return (v2 - v1).norm();
}

double triangleArea(
	const Eigen::Vector3d & v1, 
	const Eigen::Vector3d & v2, 
	const Eigen::Vector3d & v3)
{
	double a = segmentLength(v1, v2);
	double b = segmentLength(v2, v3);
	double c = segmentLength(v3, v1);

	double p = 0.5 * (a + b + c);

	return sqrt(p * (p - a) * (p - b) * (p - c));
}

Eigen::Vector3d triangleNormal(const Eigen::Vector3d & v1, const Eigen::Vector3d & v2, const Eigen::Vector3d & v3)
{
	return ((v2 - v1).cross(v3 - v1)).normalized();
}
