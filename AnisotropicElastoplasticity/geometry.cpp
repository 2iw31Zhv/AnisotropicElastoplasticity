#include "geometry.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

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

void GramSchmidtOrthonomalization(Eigen::Matrix3d & Q, 
	Eigen::Matrix3d & R, 
	const Eigen::Matrix3d & matrix)
{
	const Vector3d& d1 = matrix.col(0);
	const Vector3d& d2 = matrix.col(1);
	const Vector3d& d3 = matrix.col(2);

	Vector3d q1 = d1.normalized();
	double r11 = d1.norm();

	double r12 = d2.dot(q1);
	Vector3d q2 = d2 - r12 * q1;
	double r22 = q2.norm();
	q2.normalize();

	double r13 = d3.dot(q1);
	double r23 = d3.dot(q2);

	Vector3d q3 = d3 - r13 * q1 - r23 * q2;
	double r33 = q3.norm();
	q3.normalize();


	Q.col(0) = q1;
	Q.col(1) = q2;
	Q.col(2) = q3;

	R << r11, r12, r13,
		0.0, r22, r23,
		0.0, 0.0, r33;
}
