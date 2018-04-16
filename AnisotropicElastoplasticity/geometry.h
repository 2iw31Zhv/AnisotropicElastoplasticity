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

