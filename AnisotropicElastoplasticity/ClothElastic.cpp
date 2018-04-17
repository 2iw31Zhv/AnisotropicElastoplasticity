#include "ClothElastic.h"
#include "TriangleCondition.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>
#include <iostream>

using namespace std;
using namespace Eigen;


#define CONSIDER_STRETCH
//#define CONSIDER_SHEAR
//#define CONSIDER_BEND


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
	std::vector< Eigen::Triplet<double> > * dampingCoeffToVelocities
)
{
	double elasticEnergy = 0.0;

	if (elasticForce != nullptr)
	{
		elasticForce->resize(V.rows(), 3);
		elasticForce->setZero();
	}

	if (dampingForce != nullptr)
	{
		dampingForce->resize(V.rows(), 3);
		dampingForce->setZero();
	}

	for (int i = 0; i < F.rows(); ++i)
	{
		const Vector3d& v1 = V.row(F(i, 0));
		const Vector3d& v2 = V.row(F(i, 1));
		const Vector3d& v3 = V.row(F(i, 2));
		const Matrix2d& inverseM = inverseMetrics[i];
		double area = areas[i];
		double bu = 1.0;
		double bv = 1.0;

#ifdef CONSIDER_STRETCH
		Vector2d stretchCondition;
		computeStretchCondition(area, v1, v2, v3, inverseM, bu, bv, stretchCondition);
		double stretchEnergy = 0.5 * stretchStiffness * stretchCondition.squaredNorm();
		elasticEnergy += stretchEnergy;
#endif

#ifdef CONSIDER_SHEAR
		double shearCondition;
		computeShearCondition(area, v1, v2, v3, inverseM, shearCondition);
		double shearEnergy = 0.5 * shearStiffness * shearCondition * shearCondition;
		elasticEnergy += shearEnergy;
#endif

#ifdef CONSIDER_BEND

		auto ADD_BEND_ENERGY = [&](const Vector3d& vi, const Vector3d& vj, const Vector3d& vk,
			const Vector3d& dual_vi) {
			double bendCondition;
			computeBendCondition(vi, vj, vk, dual_vi, bendCondition);

			double bendEnergy = 0.5 * bendStiffness * bendCondition * bendCondition;
			elasticEnergy += bendEnergy;
		};

		if (faceWing(i, 0) != -1)
		{
			ADD_BEND_ENERGY(v1, v2, v3, V.row(faceWing(i, 0)));
		}

		if (faceWing(i, 1) != -1)
		{
			ADD_BEND_ENERGY(v2, v3, v1, V.row(faceWing(i, 1)));
		}

		if (faceWing(i, 2) != -1)
		{
			ADD_BEND_ENERGY(v3, v1, v2, V.row(faceWing(i, 2)));
		}
#endif

#ifdef CONSIDER_STRETCH
		if (elasticForce != nullptr)
		{
			Matrix<double, 2, 9> stretchConditionJacobian;
			computeStretchConditionJacobian(area, v1, v2, v3, inverseM,
				bu, bv, stretchConditionJacobian);
			
			elasticForce->row(F(i, 0)) += -stretchStiffness
				* stretchCondition.transpose() * stretchConditionJacobian.block(0, 0, 2, 3);
			elasticForce->row(F(i, 1)) += -stretchStiffness
				* stretchCondition.transpose() * stretchConditionJacobian.block(0, 3, 2, 3);
			elasticForce->row(F(i, 2)) += -stretchStiffness
				* stretchCondition.transpose() * stretchConditionJacobian.block(0, 6, 2, 3);

			Vector2d stretchConditionVelocity =
				stretchConditionJacobian.block(0, 0, 2, 3)
				* velocities.row(F(i, 0)).transpose()
				+ stretchConditionJacobian.block(0, 3, 2, 3)
				* velocities.row(F(i, 1)).transpose()
				+ stretchConditionJacobian.block(0, 6, 2, 3)
				* velocities.row(F(i, 2)).transpose();

			if (dampingForce != nullptr)
			{
				dampingForce->row(F(i, 0)) += -stretchStiffness * dampingStiffratio
					* stretchConditionVelocity.transpose()
					* stretchConditionJacobian.block(0, 0, 2, 3);
				dampingForce->row(F(i, 1)) += -stretchStiffness * dampingStiffratio
					* stretchConditionVelocity.transpose()
					* stretchConditionJacobian.block(0, 3, 2, 3);
				dampingForce->row(F(i, 2)) += -stretchStiffness * dampingStiffratio
					* stretchConditionVelocity.transpose()
					* stretchConditionJacobian.block(0, 6, 2, 3);
			}

			if (elasticCoeff != nullptr)
			{
				Matrix<double, 9, 9> stretchConditionHessian_Row1;
				Matrix<double, 9, 9> stretchConditionHessian_Row2;

				computeStretchConditionHessian(area, v1, v2, v3, inverseM, bu, bv,
					stretchConditionHessian_Row1,
					stretchConditionHessian_Row2);

				Matrix<double, 9, 9> K;
				K = -stretchStiffness * ((stretchConditionJacobian.transpose()
					* stretchConditionJacobian)
					+ stretchConditionHessian_Row1 * stretchCondition[0]
					+ stretchConditionHessian_Row2 * stretchCondition[1]);

				Matrix<double, 9, 9> Kdx;
				Kdx = -stretchStiffness * dampingStiffratio * (
					stretchConditionHessian_Row1 * stretchConditionVelocity[0]
					+ stretchConditionHessian_Row2 * stretchConditionVelocity[1]);
				

				Matrix<double, 9, 9> Kdv;
				Kdv = -stretchStiffness * dampingStiffratio * (
					stretchConditionJacobian.transpose()
					* stretchConditionJacobian);
				
				for (int ii = 0; ii < 3; ++ii)
				{
					for (int jj = 0; jj < 3; ++jj)
					{
						const Matrix<double, 3, 3>& blockK = K.block(3 * ii, 3 * jj, 3, 3);
						const Matrix<double, 3, 3>& blockKdx = Kdx.block(3 * ii, 3 * jj, 3, 3);
						const Matrix<double, 3, 3>& blockKdv = Kdv.block(3 * ii, 3 * jj, 3, 3);

						int idRowBase = F(i, ii);
						int idColBase = F(i, jj);
						
						for (int iii = 0; iii < 3; ++iii)
						{
							for (int jjj = 0; jjj < 3; ++jjj)
							{
								elasticCoeff->push_back(Triplet<double>(3 * idRowBase + iii,
									3 * idColBase + jjj,
									blockK(iii, jjj)));
								dampingCoeffToPositions->push_back(Triplet<double>(3 * idRowBase + iii,
									3 * idColBase + jjj,
									blockKdx(iii, jjj)));
								dampingCoeffToVelocities->push_back(Triplet<double>(3 * idRowBase + iii,
									3 * idColBase + jjj,
									blockKdv(iii, jjj)));
							}
						}
					}
				}

			}
#endif

#ifdef CONSIDER_SHEAR
			Matrix<double, 1, 9> shearConditionJacobian;
			computeShearConditionJacobian(area, v1, v2, v3, inverseM, shearConditionJacobian);

			elasticForce->row(F(i, 0)) += -shearStiffness
				* shearCondition * shearConditionJacobian.block(0, 0, 1, 3);
			elasticForce->row(F(i, 1)) += -shearStiffness
				* shearCondition * shearConditionJacobian.block(0, 3, 1, 3);
			elasticForce->row(F(i, 2)) += -shearStiffness
				* shearCondition * shearConditionJacobian.block(0, 6, 1, 3);

			double shearConditionVelocity;
			shearConditionVelocity
				= shearConditionJacobian.block(0, 0, 1, 3).dot(velocities.row(F(i, 0)))
				+ shearConditionJacobian.block(0, 3, 1, 3).dot(velocities.row(F(i, 1)))
				+ shearConditionJacobian.block(0, 6, 1, 3).dot(velocities.row(F(i, 2)));

			if (dampingForce != nullptr)
			{
				dampingForce->row(F(i, 0)) += -shearStiffness * dampingStiffratio
					* shearConditionJacobian.block(0, 0, 1, 3) * shearConditionVelocity;
				dampingForce->row(F(i, 1)) += -shearStiffness * dampingStiffratio
					* shearConditionJacobian.block(0, 3, 1, 3) * shearConditionVelocity;
				dampingForce->row(F(i, 2)) += -shearStiffness * dampingStiffratio
					* shearConditionJacobian.block(0, 6, 1, 3) * shearConditionVelocity;
			}

			if (elasticCoeff != nullptr)
			{
				Matrix<double, 9, 9> shearConditionHessian;
				computeShearConditionHessian(area, v1, v2, v3, inverseM, shearConditionHessian);

				Matrix<double, 9, 9> K;
				K = -shearStiffness * ((shearConditionJacobian.transpose()
					* shearConditionJacobian)
					+ shearConditionHessian * shearCondition);

				Matrix<double, 9, 9> Kdv;
				Matrix<double, 9, 9> Kdx;

				Kdx = -shearStiffness * dampingStiffratio
					* shearConditionHessian * shearConditionVelocity;
				Kdv = -shearStiffness * dampingStiffratio
					* shearConditionJacobian.transpose()
					* shearConditionJacobian;

				for (int ii = 0; ii < 3; ++ii)
				{
					for (int jj = 0; jj < 3; ++jj)
					{
						const Matrix<double, 3, 3>& blockK = K.block(3 * ii, 3 * jj, 3, 3);
						const Matrix<double, 3, 3>& blockKdv = Kdv.block(3 * ii, 3 * jj, 3, 3);
						const Matrix<double, 3, 3>& blockKdx = Kdx.block(3 * ii, 3 * jj, 3, 3);

						int idRowBase = F(i, ii);
						int idColBase = F(i, jj);

						for (int iii = 0; iii < 3; ++iii)
						{
							for (int jjj = 0; jjj < 3; ++jjj)
							{
								elasticCoeff->push_back(Triplet<double>(3 * idRowBase + iii,
									3 * idColBase + jjj,
									blockK(iii, jjj)));
								dampingCoeffToPositions->push_back(Triplet<double>(3 * idRowBase + iii,
									3 * idColBase + jjj,
									blockKdx(iii, jjj)));
								dampingCoeffToVelocities->push_back(Triplet<double>(3 * idRowBase + iii,
									3 * idColBase + jjj,
									blockKdv(iii, jjj)));
							}
						}
					}
				}

			}
#endif

#ifdef CONSIDER_BEND
			auto ADD_BEND_FORCE = [&](const Vector3d& vi, const Vector3d& vj, const Vector3d& vk,
				const Vector3d& dual_vi, int id_i , int id_j, int id_k, int id_di) {

				double bendCondition;
				computeBendCondition(vi, vj, vk, dual_vi, bendCondition);
				if (isnan(bendCondition))
				{
					cerr << "bendCondition is nan!\n";
					exit(1);
				}

				Matrix<double, 1, 12> bendConditionJacobian;
				computeBendConditionJacobian(vi, vj, vk, dual_vi, bendConditionJacobian);

				if (bendConditionJacobian.hasNaN())
				{
					cerr << "bendConditionJacobian.hasNaN()!\n";
					exit(1);
				}
				elasticForce->row(id_i) += - bendStiffness
					* bendCondition * bendConditionJacobian.block(0, 0, 1, 3);
				elasticForce->row(id_j) += -bendStiffness
					* bendCondition * bendConditionJacobian.block(0, 3, 1, 3);
				elasticForce->row(id_k) += -bendStiffness
					* bendCondition * bendConditionJacobian.block(0, 6, 1, 3);
				elasticForce->row(id_di) += -bendStiffness
					* bendCondition * bendConditionJacobian.block(0, 9, 1, 3);
			};

			if (faceWing(i, 0) != -1)
			{
				ADD_BEND_FORCE(v1, v2, v3, V.row(faceWing(i, 0)), 
					F(i, 0), 
					F(i, 1), 
					F(i, 2), 
					faceWing(i, 0));
			}

			if (faceWing(i, 1) != -1)
			{
				ADD_BEND_FORCE(v2, v3, v1, V.row(faceWing(i, 1)), 
					F(i, 1), 
					F(i, 2), 
					F(i, 0), 
					faceWing(i, 1));
			}

if (faceWing(i, 2) != -1)
{
	ADD_BEND_FORCE(v3, v1, v2, V.row(faceWing(i, 2)),
		F(i, 2),
		F(i, 0),
		F(i, 1),
		faceWing(i, 2));
}

if (elasticCoeff != nullptr)
{
	auto ADD_BEND_COEFF = [&](const Vector3d& vi, const Vector3d& vj, const Vector3d& vk,
		const Vector3d& dual_vi, int id_i, int id_j, int id_k, int id_di)
	{
		double bendCondition;
		computeBendCondition(vi, vj, vk, dual_vi, bendCondition);
		Matrix<double, 1, 12> bendConditionJacobian;
		computeBendConditionJacobian(vi, vj, vk, dual_vi, bendConditionJacobian);
		Matrix<double, 12, 12> bendConditionHessian;
		computeBendConditionHessian(vi, vj, vk, dual_vi, bendConditionHessian);

		Matrix<double, 12, 12> K;
		K = - bendStiffness * ((bendConditionJacobian.transpose()
			* bendConditionJacobian)
			+ bendConditionHessian * bendCondition
			);

		if (K.hasNaN())
		{
			cerr << "K.hasNaN()!\n";
			exit(1);
		}

		int ids[4] = { id_i, id_j, id_k, id_di };

		for (int ii = 0; ii < 4; ++ii)
		{
			for (int jj = 0; jj < 4; ++jj)
			{
				int idRowBase = ids[ii];
				int idColBase = ids[jj];

				const Matrix3d& blockK = K.block(3 * ii, 3 * jj, 3, 3);

				for (int iii = 0; iii < 3; ++iii)
				{
					for (int jjj = 0; jjj < 3; ++jjj)
					{
						elasticCoeff->push_back(Triplet<double>(3 * idRowBase + iii,
							3 * idColBase + jjj,
							blockK(iii, jjj)));
					}
				}
			}
		}

	};

	if (faceWing(i, 0) != -1)
	{
		ADD_BEND_COEFF(v1, v2, v3,
			V.row(faceWing(i, 0)),
			F(i, 0),
			F(i, 1),
			F(i, 2),
			faceWing(i, 0));
	}
	if (faceWing(i, 1) != -1)
	{
		ADD_BEND_COEFF(v2, v3, v1,
			V.row(faceWing(i, 1)),
			F(i, 1),
			F(i, 2),
			F(i, 0),
			faceWing(i, 1));
	}
	if (faceWing(i, 2) != -1)
	{
		ADD_BEND_COEFF(v3, v1, v2,
			V.row(faceWing(i, 2)),
			F(i, 2),
			F(i, 0),
			F(i, 1),
			faceWing(i, 2));
	}
}

#endif			
		}
	}

	return elasticEnergy;
}

bool hasDrasticMetricChange(
	const Eigen::MatrixX3d & oldV,
	const Eigen::MatrixX3d & V,
	const Eigen::MatrixX3i& F,
	const std::vector<Eigen::Matrix2d>& inverseMetrics)
{
	for (int i = 0; i < F.rows(); ++i)
	{
		const Vector3d& v1 = V.row(F(i, 0));
		const Vector3d& v2 = V.row(F(i, 1));
		const Vector3d& v3 = V.row(F(i, 2));
		if (v1.hasNaN() || v2.hasNaN() || v3.hasNaN())
		{
			return true;
		}
		const Vector3d& oldv1 = oldV.row(F(i, 0));
		const Vector3d& oldv2 = oldV.row(F(i, 1));
		const Vector3d& oldv3 = oldV.row(F(i, 2));
		const Matrix2d& inverseM = inverseMetrics[i];

		Vector3d wu, wv, oldWu, oldWv;

		triangleMetric(v1, v2, v3, inverseM, wu, wv);
		triangleMetric(v1, v2, v3, inverseM, oldWu, oldWv);



		if ((wu - oldWu).squaredNorm() > 1e-3 ||
			(wv - oldWv).squaredNorm() > 1e-3)
		{
			cerr << (wu - oldWu).squaredNorm() << ", "
				<< (wv - oldWv).squaredNorm() << endl;
			return true;
		}
	}

	return false;
}
