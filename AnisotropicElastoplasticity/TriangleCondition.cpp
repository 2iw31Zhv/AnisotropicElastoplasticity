#include "TriangleCondition.h"
#include <Eigen/Dense>
#include "LinearAlgebra.h"
#include "MathExtensions.h"

#include <iostream>

using namespace std;

using namespace Eigen;

void triangleMetric(const Eigen::Vector3d & xi, 
	const Eigen::Vector3d & xj, 
	const Eigen::Vector3d & xk, 
	const Eigen::Matrix2d& inverseM,
	Eigen::Vector3d & wu, 
	Eigen::Vector3d & wv)
{
	Vector3d dx1 = xj - xi;
	Vector3d dx2 = xk - xi;

	wu = inverseM(0, 0) * dx1 + inverseM(1, 0) * dx2;
	wv = inverseM(0, 1) * dx1 + inverseM(1, 1) * dx2;
}

void triangleMetricJacobian_(const Eigen::Vector3d & xi,
	const Eigen::Vector3d & xj,
	const Eigen::Vector3d & xk,
	const Eigen::Matrix2d& inverseM,
	Eigen::Matrix<double, 3, 9>& wuJacobian,
	Eigen::Matrix<double, 3, 9>& wvJacobian)
{    
	// wu = m00 * (xj - xi) + m10 * (xk - xi)
	// wv = m01 * (xj - xi) + m11 * (xk - xi)

	wuJacobian.block(0, 0, 3, 3) = -(inverseM(0, 0) + inverseM(1, 0)) * Matrix3d::Identity();
	wuJacobian.block(0, 3, 3, 3) = inverseM(0, 0) * Matrix3d::Identity();
	wuJacobian.block(0, 6, 3, 3) = inverseM(1, 0) * Matrix3d::Identity();

	wvJacobian.block(0, 0, 3, 3) = -(inverseM(0, 1) + inverseM(1, 1)) * Matrix3d::Identity();
	wvJacobian.block(0, 3, 3, 3) = inverseM(0, 1) * Matrix3d::Identity();
	wvJacobian.block(0, 6, 3, 3) = inverseM(1, 1) * Matrix3d::Identity();
}

void computeStretchCondition_(double area, 
	const Eigen::Vector3d & wu,
	const Eigen::Vector3d & wv, 
	double bu, double bv, 
	Eigen::Vector2d & stretchCondition)
{
	stretchCondition(0) = area * (wu.norm() - bu);
	stretchCondition(1) = area * (wv.norm() - bv);
}

void computeShearCondition_(double area, 
	const Eigen::Vector3d & wu, 
	const Eigen::Vector3d & wv, 
	double & shearCondition)
{
	shearCondition = area * wu.dot(wv);
}

void computeStretchCondition(double area, 
	const Eigen::Vector3d & xi, 
	const Eigen::Vector3d & xj, 
	const Eigen::Vector3d & xk, 
	const Eigen::Matrix2d& inverseM,
	double bu, 
	double bv, 
	Eigen::Vector2d & stretchCondition)
{
	Vector3d wu, wv;
	triangleMetric(xi, xj, xk, inverseM, wu, wv);
	computeStretchCondition_(area, wu, wv, bu, bv, stretchCondition);
}

void computeStretchConditionJacobian(double area, 
	const Eigen::Vector3d & xi, 
	const Eigen::Vector3d & xj, 
	const Eigen::Vector3d & xk, 
	const Eigen::Matrix2d& inverseM,
	double bu, 
	double bv, 
	Eigen::Matrix<double, 2, 9>& stretchConditionJacobian)
{
	Vector3d wu, wv;
	triangleMetric(xi, xj, xk, inverseM, wu, wv);
	Matrix<double, 3, 9> wuJacobian, wvJacobian;
	triangleMetricJacobian_(xi, xj, xk, inverseM, wuJacobian, wvJacobian);

	stretchConditionJacobian.row(0) = area / wu.norm() * wu.transpose() * wuJacobian;
	stretchConditionJacobian.row(1) = area / wv.norm() * wv.transpose() * wvJacobian;
}

void computeStretchConditionHessian(double area,
	const Eigen::Vector3d & xi,
	const Eigen::Vector3d & xj,
	const Eigen::Vector3d & xk,
	const Eigen::Matrix2d & inverseM,
	double bu, double bv,
	Eigen::Matrix<double, 9, 9>& stretchConditionHessian_Row1,
	Eigen::Matrix<double, 9, 9>& stretchConditionHessian_Row2)
{
	Vector3d wu, wv;
	triangleMetric(xi, xj, xk, inverseM, wu, wv);
	Matrix<double, 3, 9> wuJacobian, wvJacobian;
	triangleMetricJacobian_(xi, xj, xk, inverseM, wuJacobian, wvJacobian);

	stretchConditionHessian_Row1 =
		-area * wuJacobian.transpose() * wu * wu.transpose() / pow(wu.norm(), 3)
		* wuJacobian;
	stretchConditionHessian_Row2 =
		-area * wvJacobian.transpose() * wv * wv.transpose() / pow(wv.norm(), 3)
		* wvJacobian;
}

void computeShearCondition(
	double area, 
	const Eigen::Vector3d & xi, 
	const Eigen::Vector3d & xj, 
	const Eigen::Vector3d & xk, 
	const Eigen::Matrix2d& inverseM,
	double & shearCondition)
{
	Vector3d wu, wv;
	triangleMetric(xi, xj, xk, inverseM, wu, wv);
	computeShearCondition_(area, wu, wv, shearCondition);
}

void computeShearConditionJacobian(
	double area, 
	const Eigen::Vector3d & xi, 
	const Eigen::Vector3d & xj, 
	const Eigen::Vector3d & xk, 
	const Eigen::Matrix2d& inverseM,
	Eigen::Matrix<double, 1, 9>& shearConditionJacobian)
{
	Vector3d wu, wv;
	triangleMetric(xi, xj, xk, inverseM, wu, wv);
	Matrix<double, 3, 9> wuJacobian, wvJacobian;
	triangleMetricJacobian_(xi, xj, xk, inverseM, wuJacobian, wvJacobian);

	shearConditionJacobian = area * (wv.transpose() * wuJacobian + wu.transpose() * wvJacobian);
}

void computeShearConditionHessian(double area, 
	const Eigen::Vector3d & xi, 
	const Eigen::Vector3d & xj, 
	const Eigen::Vector3d & xk, 
	const Eigen::Matrix2d & inverseM, 
	Eigen::Matrix<double, 9, 9>& shearConditionHessian)
{
	Vector3d wu, wv;
	triangleMetric(xi, xj, xk, inverseM, wu, wv);
	Matrix<double, 3, 9> wuJacobian, wvJacobian;
	triangleMetricJacobian_(xi, xj, xk, inverseM, wuJacobian, wvJacobian);

	shearConditionHessian = area * (wvJacobian.transpose() * wuJacobian + wuJacobian.transpose() * wvJacobian);
}

void computeBendCondition(const Eigen::Vector3d & xi, 
	const Eigen::Vector3d & xj, 
	const Eigen::Vector3d & xk, 
	const Eigen::Vector3d & dual_xi, 
	double & bendCondition)
{
	Vector3d n1 = (xj - xi).cross(xk - xi); n1.normalize();
	Vector3d n2 = (xk - dual_xi).cross(xj - dual_xi); n2.normalize();

	double x = n1.dot(n2);

	if (x < -1.0)
	{
		bendCondition = MathExtensions::Pi;
	}
	else if (1.0 < x)
	{
		bendCondition = 0.0;
	}
	else
	{
		bendCondition = MathExtensions::acos(x);
	}

	if (isnan(bendCondition))
	{
		cerr << "xi: " << xi.x() << ", " << xi.y() << ", " << xi.z() << endl;
		cerr << "xj: " << xj.x() << ", " << xj.y() << ", " << xj.z() << endl;
		cerr << "xk: " << xk.x() << ", " << xk.y() << ", " << xk.z() << endl;
		cerr << "dual_xi: " << dual_xi.x() << ", " << dual_xi.y() << ", " << dual_xi.z() << endl;
	}
}

void computeFacetNormalJacobian_(
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	Eigen::Matrix<double, 3, 9>& normalJacobian)
{
	Vector3d ejk = xk - xj;
	Vector3d eki = xi - xk;
	Vector3d eij = xj - xi;
	normalJacobian.block(0, 0, 3, 3) = crossMatrix(ejk);
	normalJacobian.block(0, 3, 3, 3) = crossMatrix(eki);
	normalJacobian.block(0, 6, 3, 3) = crossMatrix(eij);
}

void computeInnerProductJacobian_(
	const Eigen::Vector3d & xi,
	const Eigen::Vector3d & xj,
	const Eigen::Vector3d & xk,
	const Eigen::Vector3d & dual_xi,
	Eigen::Matrix<double, 1, 12>& innerProductJacobian)
{
	Vector3d n1 = (xj - xi).cross(xk - xi); n1.normalize();
	Vector3d n2 = (xk - dual_xi).cross(xj - dual_xi); n2.normalize();
	double innerProduct = n1.dot(n2);

	Matrix<double, 3, 9> n1Jacobian, n2Jacobian;
	computeFacetNormalJacobian_(xi, xj, xk, n1Jacobian);
	computeFacetNormalJacobian_(dual_xi, xk, xj, n2Jacobian);


	Matrix<double, 3, 12> dn1, dn2;
	dn1.setZero();
	dn2.setZero();
	dn1.block(0, 0, 3, 9) = n1Jacobian;
	dn2.block(0, 3, 3, 3) = n2Jacobian.block(0, 6, 3, 3);
	dn2.block(0, 6, 3, 3) = n2Jacobian.block(0, 3, 3, 3);
	dn2.block(0, 9, 3, 3) = n2Jacobian.block(0, 0, 3, 3);

	innerProductJacobian = n1.transpose() * dn2 + n2.transpose() * dn1;
}

void computeCrossProductJacobian_(
	const Eigen::Vector3d & xi,
	const Eigen::Vector3d & xj,
	const Eigen::Vector3d & xk,
	const Eigen::Vector3d & dual_xi,
	Eigen::Matrix<double, 3, 12>& crossProductJacobian)
{
	Vector3d n1 = (xj - xi).cross(xk - xi); n1.normalize();
	Vector3d n2 = (xk - dual_xi).cross(xj - dual_xi); n2.normalize();

	Matrix<double, 3, 3> crossProductToN1 = -crossMatrix(n2);
	Matrix<double, 3, 3> crossProductToN2 = crossMatrix(n1);
	
	Matrix<double, 3, 9> n1Jacobian, n2Jacobian;
	computeFacetNormalJacobian_(xi, xj, xk, n1Jacobian);
	computeFacetNormalJacobian_(dual_xi, xk, xj, n2Jacobian);


	Matrix<double, 3, 12> dn1, dn2;
	dn1.setZero();
	dn2.setZero();
	dn1.block(0, 0, 3, 9) = n1Jacobian;
	dn2.block(0, 3, 3, 3) = n2Jacobian.block(0, 6, 3, 3);
	dn2.block(0, 6, 3, 3) = n2Jacobian.block(0, 3, 3, 3);
	dn2.block(0, 9, 3, 3) = n2Jacobian.block(0, 0, 3, 3);
	
	crossProductJacobian = crossProductToN1 * dn1
		+ crossProductToN2 * dn2;
}

void computeBendConditionJacobian(
	const Eigen::Vector3d & xi, 
	const Eigen::Vector3d & xj, 
	const Eigen::Vector3d & xk, 
	const Eigen::Vector3d & dual_xi, 
	Eigen::Matrix<double, 1, 12>& bendConditionJacobian)
{
	Vector3d n1 = (xj - xi).cross(xk - xi); n1.normalize();
	Vector3d n2 = (xk - dual_xi).cross(xj - dual_xi); n2.normalize();
	double innerProduct = n1.dot(n2);

	Vector3d crossProduct = n1.cross(n2);
	Vector3d edge = xj - xk;

	edge.normalize();
	Eigen::Matrix<double, 1, 12> innerProductJacobian;
	computeInnerProductJacobian_(xi, xj, xk, dual_xi, innerProductJacobian);
	Eigen::Matrix<double, 3, 12> crossProductJacobian;
	computeCrossProductJacobian_(xi, xj, xk, dual_xi, crossProductJacobian);

	if (fabs(1 - innerProduct * innerProduct) >= 1e-5)
	{
		bendConditionJacobian = MathExtensions::dacos(innerProduct) * innerProductJacobian;
	}
	else
	{
		bendConditionJacobian = MathExtensions::dasin(crossProduct.norm())
			* edge.transpose()
			* crossProductJacobian;  
	}
	
	if (bendConditionJacobian.hasNaN())
	{
		cerr << "bendConditionJacobian.hasNaN()!\n";
		exit(1);
	}
}

void computeFacetNormalHessian_(
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	Eigen::Matrix<double, 9, 9>& normalHessian_x,
	Eigen::Matrix<double, 9, 9>& normalHessian_y,
	Eigen::Matrix<double, 9, 9>& normalHessian_z)
{
	normalHessian_x.setZero();
	normalHessian_y.setZero();
	normalHessian_z.setZero();

	const int x1 = 0;
	const int y1 = 1;
	const int z1 = 2;
	const int x2 = 3;
	const int y2 = 4;
	const int z2 = 5;
	const int x3 = 6;
	const int y3 = 7;
	const int z3 = 8;

	// set x
	normalHessian_x(1, z2) = 1.0;
	normalHessian_x(1, z3) = -1.0;
	normalHessian_x(2, y2) = -1.0;
	normalHessian_x(2, y3) = 1.0;

	normalHessian_x(4, z3) = 1.0;
	normalHessian_x(4, z1) = -1.0;
	normalHessian_x(5, y3) = -1.0;
	normalHessian_x(5, y1) = 1.0;

	normalHessian_x(7, z1) = 1.0;
	normalHessian_x(7, z2) = -1.0;
	normalHessian_x(8, y1) = -1.0;
	normalHessian_x(8, y2) = 1.0;

	// set y
	normalHessian_y(0, z2) = -1.0;
	normalHessian_y(0, z3) = 1.0;
	normalHessian_y(2, x2) = 1.0;
	normalHessian_y(2, x3) = -1.0;
	
	normalHessian_y(3, z3) = -1.0;
	normalHessian_y(3, z1) = 1.0;
	normalHessian_y(5, x3) = 1.0;
	normalHessian_y(5, x1) = -1.0;

	normalHessian_y(6, z1) = -1.0;
	normalHessian_y(6, z2) = 1.0;
	normalHessian_y(8, x1) = 1.0;
	normalHessian_y(8, x2) = -1.0;

	// set z
	normalHessian_z(0, y2) = 1.0;
	normalHessian_z(0, y3) = -1.0;
	normalHessian_z(1, x2) = -1.0;
	normalHessian_z(1, x3) = 1.0;

	normalHessian_z(3, y3) = 1.0;
	normalHessian_z(3, y1) = -1.0;
	normalHessian_z(4, x3) = -1.0;
	normalHessian_z(4, x1) = 1.0;

	normalHessian_z(6, y1) = 1.0;
	normalHessian_z(6, y2) = -1.0;
	normalHessian_z(6, x1) = -1.0;
	normalHessian_z(6, x2) = 1.0;
}

void computeInnerProductHessian_(
	const Eigen::Vector3d& xi,
	const Eigen::Vector3d& xj,
	const Eigen::Vector3d& xk,
	const Eigen::Vector3d& dual_xi,
	Eigen::Matrix<double, 12, 12>& innerProductHessian
)
{
	Vector3d n1 = (xj - xi).cross(xk - xi); n1.normalize();
	Vector3d n2 = (xk - dual_xi).cross(xj - dual_xi); n2.normalize();
	double innerProduct = n1.dot(n2);

	Matrix<double, 3, 9> n1Jacobian, n2Jacobian;
	computeFacetNormalJacobian_(xi, xj, xk, n1Jacobian);
	computeFacetNormalJacobian_(dual_xi, xk, xj, n2Jacobian);


	Matrix<double, 3, 12> dn1, dn2;
	dn1.setZero();
	dn2.setZero();
	dn1.block(0, 0, 3, 9) = n1Jacobian;
	dn2.block(0, 3, 3, 3) = n2Jacobian.block(0, 6, 3, 3);
	dn2.block(0, 6, 3, 3) = n2Jacobian.block(0, 3, 3, 3);
	dn2.block(0, 9, 3, 3) = n2Jacobian.block(0, 0, 3, 3);

	Matrix<double, 9, 9> n1Hessian_x, n1Hessian_y, n1Hessian_z;
	Matrix<double, 9, 9> n2Hessian_x, n2Hessian_y, n2Hessian_z;
	computeFacetNormalHessian_(xi, xj, xk, n1Hessian_x, n1Hessian_y, n1Hessian_z);
	computeFacetNormalHessian_(dual_xi, xk, xj, n2Hessian_x, n2Hessian_y, n2Hessian_z);

	Matrix<double, 12, 12> ddn1_x, ddn1_y, ddn1_z;
	Matrix<double, 12, 12> ddn2_x, ddn2_y, ddn2_z;
	ddn1_x.setZero(); ddn1_y.setZero(); ddn1_z.setZero();
	ddn1_x.block(0, 0, 9, 9) = n1Hessian_x;
	ddn1_y.block(0, 0, 9, 9) = n1Hessian_y;
	ddn1_z.block(0, 0, 9, 9) = n1Hessian_z;

	auto SET_MATRIX = [&](Matrix<double, 12, 12>& ddn2,
		const Matrix<double, 9, 9>& n2Hessian)
	{
		ddn2.setZero();
		ddn2.block(3, 3, 3, 3) = n2Hessian.block(6, 6, 3, 3);
		ddn2.block(3, 6, 3, 3) = n2Hessian.block(6, 3, 3, 3);
		ddn2.block(3, 9, 3, 3) = n2Hessian.block(6, 0, 3, 3);
		ddn2.block(6, 3, 3, 3) = n2Hessian.block(3, 6, 3, 3);
		ddn2.block(6, 6, 3, 3) = n2Hessian.block(3, 3, 3, 3);
		ddn2.block(6, 9, 3, 3) = n2Hessian.block(3, 0, 3, 3);
		ddn2.block(9, 3, 3, 3) = n2Hessian.block(0, 6, 3, 3);
		ddn2.block(9, 6, 3, 3) = n2Hessian.block(0, 3, 3, 3);
		ddn2.block(9, 9, 3, 3) = n2Hessian.block(0, 0, 3, 3);
	};
	SET_MATRIX(ddn2_x, n2Hessian_x);
	SET_MATRIX(ddn2_y, n2Hessian_y);
	SET_MATRIX(ddn2_z, n2Hessian_z);

	innerProductHessian = 2 * dn1.transpose() * dn2
		+ n1.x() * ddn1_x + n1.y() * ddn1_y + n1.z() * ddn1_z
		+ n2.x() * ddn2_x + n2.y() * ddn2_y + n2.z() * ddn2_z;
}

void computeBendConditionHessian(
	const Eigen::Vector3d & xi, 
	const Eigen::Vector3d & xj, 
	const Eigen::Vector3d & xk, 
	const Eigen::Vector3d & dual_xi, 
	Eigen::Matrix<double, 12, 12>& bendConditionHessian)
{
	Vector3d n1 = (xj - xi).cross(xk - xi); n1.normalize();
	Vector3d n2 = (xk - dual_xi).cross(xj - dual_xi); n2.normalize();
	double innerProduct = n1.dot(n2);
	Vector3d crossProduct = n1.cross(n2);

	Matrix<double, 1, 12> innerProductJacobian;
	computeInnerProductJacobian_(xi, xj, xk, dual_xi, innerProductJacobian);

	Matrix<double, 12, 12> innerProductHessian;
	computeInnerProductHessian_(xi, xj, xk, dual_xi, innerProductHessian);

	if (fabs(1 - innerProduct * innerProduct) >= 1e-5)
	{
		bendConditionHessian = MathExtensions::dacos(innerProduct) * innerProductHessian
			+ MathExtensions::ddacos(innerProduct) * innerProductJacobian.transpose() * innerProductJacobian;
	}
	else
	{
		double x = innerProduct;
		bendConditionHessian = -1.0 / sqrt(1.0 - x * x + 1e-5)* innerProductHessian
			 - x / pow(1.0 - x * x + 1e-5, 1.5) *  innerProductJacobian.transpose() * innerProductJacobian;
	}
	
	if (bendConditionHessian.hasNaN())
	{
		cerr << "bendConditionHessian.hasNaN()!\n";
		exit(1);

	}
}

