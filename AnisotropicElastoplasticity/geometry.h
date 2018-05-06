#pragma once

#include <Eigen/Core>

double segmentLength(const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2);

double triangleArea(const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2,
	const Eigen::Vector3d& v3);

Eigen::Vector3d triangleNormal(const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2,
	const Eigen::Vector3d& v3);

void GramSchmidtOrthonomalization(
	Eigen::Matrix3d& Q,
	Eigen::Matrix3d& R,
	const Eigen::Matrix3d& matrix);

Eigen::Matrix2d inverseR(const Eigen::Matrix2d& R);

// with BUG
void PolarDecompositionR(
	Eigen::Matrix2d& rotation,
	Eigen::Matrix2d& symmetry,
	const Eigen::Matrix2d& R);