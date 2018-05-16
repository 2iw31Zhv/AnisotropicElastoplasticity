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

#define POW2(x) ((x) * (x))
#define POW4(x) ((x) * (x) * (x) * (x))

Eigen::Matrix2d inverseR(const Eigen::Matrix2d& R)
{
	Matrix2d invR;
	invR << 1.0 / R(0, 0), -R(0, 1) / R(0, 0) / R(1, 1),
		0.0, 1.0 / R(1, 1);
	return invR;
}

void PolarDecompositionR(
	Eigen::Matrix2d& rotation,
	Eigen::Matrix2d& symmetry,
	const Eigen::Matrix2d& R)
{
	double r11 = R(0, 0),
		r12 = R(0, 1), r22 = R(1, 1);
	
	double sumSqr = POW2(r11) + POW2(r12) + POW2(r22);
	double delta = POW2(sumSqr) - 4.0 * POW2(r11) * POW2(r22);

	Vector2d Sigma;
	Sigma[0] = sqrt(0.5 * (sumSqr - sqrt(delta)));
	Sigma[1] = sqrt(0.5 * (sumSqr + sqrt(delta)));

	double pmSqr = POW2(r11) + POW2(r12) - POW2(r22);
	double PmD = pmSqr - sqrt(delta);
	double PpD = pmSqr + sqrt(delta);

	double mD = sqrt(4 * POW2(12) * POW2(22) + POW2(PmD));
	double pD = sqrt(4 * POW2(12) * POW2(22) + POW2(PpD));

	Matrix2d U;
	U << PmD / mD, PpD / pD,
		2.0 * r11 * r12 / mD, 2.0 * r11 * r12 / pD;

	rotation = U * Sigma.cwiseInverse().asDiagonal() * U.transpose() * R;
	symmetry = rotation.transpose() * R;
}
