#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

double computeClothElasticEnergy(
	const Eigen::MatrixX3d& velocities,
	const Eigen::MatrixX3d& V,
	const Eigen::MatrixX3i& F,
	const Eigen::MatrixX3i& faceWing,
	const std::vector< Eigen::Matrix2d >& inverseMetrics,
	const Eigen::VectorXd& areas,
	double stretchStiffness,
	double shearStiffness,
	double bendStiffness,
	Eigen::MatrixX3d * elasticForce,
	std::vector< Eigen::Triplet<double> > * elasticCoeff,
	double dampingStiffratio,
	Eigen::MatrixX3d * dampingForce,
	std::vector< Eigen::Triplet<double> > * dampingCoeffToPositions,
	std::vector< Eigen::Triplet<double> > * dampingCoeffToVelocities);

bool hasDrasticMetricChange(
	const Eigen::MatrixX3d& oldV,
	const Eigen::MatrixX3d& V,
	const Eigen::MatrixX3i& F,
	const std::vector< Eigen::Matrix2d >& inverseMetrics);