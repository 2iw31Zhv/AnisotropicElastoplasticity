#pragma once
#include <Eigen/Core>

void triangleMetric(const Eigen::Vector3d & xi,
	const Eigen::Vector3d & xj,
	const Eigen::Vector3d & xk,
	const Eigen::Matrix2d& inverseM,
	Eigen::Vector3d & wu,
	Eigen::Vector3d & wv);

void computeStretchCondition(
	double area,
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	const Eigen::Matrix2d& inverseM,
	double bu, double bv,
	Eigen::Vector2d& stretchCondition);

void computeStretchConditionJacobian(
	double area,
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	const Eigen::Matrix2d& inverseM,
	double bu, double bv,
	Eigen::Matrix<double, 2, 9>& stretchConditionJacobian);

void computeStretchConditionHessian(
	double area,
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	const Eigen::Matrix2d& inverseM,
	double bu, double bv,
	Eigen::Matrix<double, 9, 9>& stretchConditionHessian_Row1,
	Eigen::Matrix<double, 9, 9>& stretchConditionHessian_Row2);

void computeShearCondition(
	double area,
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	const Eigen::Matrix2d& inverseM,
	double& shearCondition);

void computeShearConditionJacobian(
	double area,
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	const Eigen::Matrix2d& inverseM,
	Eigen::Matrix<double, 1, 9>& shearConditionJacobian);

void computeShearConditionHessian(
	double area,
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	const Eigen::Matrix2d& inverseM,
	Eigen::Matrix<double, 9, 9>& shearConditionHessian);

void computeBendCondition(
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	const Eigen::Vector3d& dual_xi,
	double& bendCondition);

void computeBendConditionJacobian(
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	const Eigen::Vector3d& dual_xi,
	Eigen::Matrix<double, 1, 12>& bendConditionJacobian);

void computeBendConditionHessian(
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	const Eigen::Vector3d& dual_xi,
	Eigen::Matrix<double, 12, 12>& bendConditionHessian);