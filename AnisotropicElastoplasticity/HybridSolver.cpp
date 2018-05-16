#include "HybridSolver.h"
#include "ParticleSystem.h"
#include "RegularGrid.h"
#include "LagrangianMesh.h"
#include <vector>
#include <igl/viewer/Viewer.h>
#include "interpolation.h"
#include <Eigen/Dense>
#include <functional>
#include "geometry.h"

using namespace std;
using namespace Eigen;
using namespace igl;

#define AFFINE_PIC

// the relaxation code, seems not to be needed

//void matrixRelaxer(
//	Eigen::MatrixX3d& matrix,
//	double tolerance = 1e-20)
//{
//	for (int p = 0; p < matrix.rows(); ++p)
//	{
//		for (int i = 0; i < 3; ++i)
//		{
//			if (fabs(matrix(p, i)) < tolerance)
//			{
//				matrix(p, i) = 0.0;
//			}
//		}
//	}
//}
//
//void matrix33Relaxer(
//	Eigen::Matrix3d& matrix,
//	double tolerance = 1e-20)
//{
//	for (int i = 0; i < 3; ++i)
//	{
//		for (int j = 0; j < 3; ++j)
//		{
//			if (i != j && fabs(matrix(i, j)) < tolerance)
//			{
//				matrix(i, j) = 0.0;
//			}
//			else if (i == j && fabs(matrix(i, j) - 1.0) < tolerance)
//			{
//				matrix(i, j) = 1.0;
//			}
//		}
//	}
//}
//
//void vector3Relaxer(
//	Eigen::Vector3d& vector,
//	double tolerance = 1e-20)
//{
//	for (int i = 0; i < 3; ++i)
//	{
//		if (fabs(vector[i]) < tolerance)
//		{
//			vector[i] = 0.0;
//		}
//	}
//}

void HybridSolver::evaluateInterpolationWeights_(
	Eigen::SparseMatrix<double>& omegas,
	Eigen::SparseMatrix<double>& domegas_1,
	Eigen::SparseMatrix<double>& domegas_2,
	Eigen::SparseMatrix<double>& domegas_3,
	const Eigen::MatrixX3d& particlePositions)
{
	vector<Triplet<double>> omegaTriplets;
	vector<Triplet<double>> domegaTriplets_1;
	vector<Triplet<double>> domegaTriplets_2;
	vector<Triplet<double>> domegaTriplets_3;

	for (int p = 0; p < particlePositions.rows(); ++p)
	{
		const Vector3d& position = particlePositions.row(p);
		const Vector3d& posRef = position - rg_->minBound();
		int i_fl = static_cast<int>(posRef[0] / rg_->h()[0]);
		int j_fl = static_cast<int>(posRef[1] / rg_->h()[1]);
		int k_fl = static_cast<int>(posRef[2] / rg_->h()[2]);

		for (int i = i_fl - 2; i <= i_fl + 2; ++i)
		{
			for (int j = j_fl - 2; j <= j_fl + 2; ++j)
			{
				for (int k = k_fl - 2; k <= k_fl + 2; ++k)
				{
					if (0 <= i && i < rg_->resolution()[0] &&
						0 <= j && j < rg_->resolution()[1] &&
						0 <= k && k < rg_->resolution()[2])
					{
						double dyadic_i = cubic_B_spline(posRef[0] / rg_->h()[0] - i);
						double dyadic_j = cubic_B_spline(posRef[1] / rg_->h()[1] - j);
						double dyadic_k = cubic_B_spline(posRef[2] / rg_->h()[2] - k);

						double ddyadic_i = Dcubic_B_spline(posRef[0] / rg_->h()[0] - i)
							/ rg_->h()[0];
						double ddyadic_j = Dcubic_B_spline(posRef[1] / rg_->h()[1] - j)
							/ rg_->h()[1];
						double ddyadic_k = Dcubic_B_spline(posRef[2] / rg_->h()[2] - k)
							/ rg_->h()[2];

						int index = rg_->toIndex(i, j, k);
						if (dyadic_i > 0 && dyadic_j > 0 && dyadic_k > 0)
						{
							omegaTriplets.push_back(Triplet<double>(p,
								index,
								dyadic_i * dyadic_j * dyadic_k));
							domegaTriplets_1.push_back(Triplet<double>(p,
								index,
								ddyadic_i * dyadic_j * dyadic_k));
							domegaTriplets_2.push_back(Triplet<double>(p,
								index,
								dyadic_i * ddyadic_j * dyadic_k));
							domegaTriplets_3.push_back(Triplet<double>(p,
								index,
								dyadic_i * dyadic_j * ddyadic_k));
						}
					}
				}
			}
		}
	}

	int Np = particlePositions.rows();
	omegas.resize(Np, rg_->masses.size());
	domegas_1.resize(Np, rg_->masses.size());
	domegas_2.resize(Np, rg_->masses.size());
	domegas_3.resize(Np, rg_->masses.size());

	omegas.setZero();
	domegas_1.setZero();
	domegas_2.setZero();
	domegas_3.setZero();

	omegas.setFromTriplets(omegaTriplets.begin(), omegaTriplets.end());
	domegas_1.setFromTriplets(domegaTriplets_1.begin(), domegaTriplets_1.end());
	domegas_2.setFromTriplets(domegaTriplets_2.begin(), domegaTriplets_2.end());
	domegas_3.setFromTriplets(domegaTriplets_3.begin(), domegaTriplets_3.end());

	//ofstream fout("omegalog.txt", ios::app);
	//VectorXd gridUnit(rg_->gridNumber());
	//gridUnit.setOnes();
	//VectorXd particleUnit(particlePositions.rows());
	//particleUnit.setOnes();

	//fout << "grid" << endl << omegas * gridUnit << endl;
	//fout << "particle" << endl << omegas.transpose() * particleUnit << endl;

}

Vector3d totalAngularMomentum(
	const Eigen::MatrixX3d& positions,
	const Eigen::MatrixX3d& momenta)
{
	Vector3d a;
	a.setZero();

	for (int p = 0; p < positions.rows(); ++p)
	{
		a += positions.row(p).cross(momenta.row(p));
	}
	return a;
}

void HybridSolver::particleToGrid_(double tolerance, bool evaluateVolumesAndDensities)
{
	rg_->masses.setZero();
	if (ps_ != nullptr)
	{
		rg_->masses += omegas_.transpose() * ps_->masses;
	}

	if (mesh_ != nullptr)
	{
		rg_->masses += vertexOmegas_.transpose() * mesh_->vertexMasses;
		rg_->masses += elementOmegas_.transpose() * mesh_->elementMasses;
	}


	MatrixX3d momenta_i;
	momenta_i.resize(rg_->gridNumber(), 3);
	momenta_i.setZero();

	if (ps_ != nullptr)
	{
		momenta_i += omegas_.transpose() * ps_->masses.asDiagonal() * ps_->velocities;
	}

	//ofstream fout("conservation_log.txt", ios::app);
	//fout << "momentum particle: ";
	//MatrixX3d vertexMomentum = mesh_->vertexMasses.asDiagonal() * mesh_->vertexVelocities;
	//MatrixX3d elementMomentum = mesh_->elementMasses.asDiagonal() * mesh_->elementVelocities;
	//Vector3d momentumParticleSum = vertexMomentum.colwise().sum()
	//	+ elementMomentum.colwise().sum();

	//fout << momentumParticleSum.transpose() << endl;

	//fout << "angular momentum particle: ";
	//fout << (totalAngularMomentum(mesh_->vertexPositions, vertexMomentum)
	//	+ totalAngularMomentum(mesh_->elementPositions, elementMomentum)).transpose() << endl;

	if (mesh_ != nullptr)
	{
		momenta_i += vertexOmegas_.transpose() * mesh_->vertexMasses.asDiagonal() * mesh_->vertexVelocities;
		momenta_i += elementOmegas_.transpose() * mesh_->elementMasses.asDiagonal() * mesh_->elementVelocities;
	}

	//fout << "momemtum grid before APIC: "
	//	<< momenta_i.colwise().sum() << endl;

	//fout << "angular mometum before APIC: "
	//	<< totalAngularMomentum(rg_->positions(), momenta_i).transpose() << endl;

#if defined(AFFINE_PIC)
	auto ADD_AFFINE_MOMENTA = [&](
		const SparseMatrix<double>& omegas,
		const VectorXd& particleMasses,
		const MatrixX3d& particleAffineMomenta_1,
		const MatrixX3d& particleAffineMomenta_2, 
		const MatrixX3d& particleAffineMomenta_3,
		const MatrixX3d& particlePositions)
	{
		

		/* the original code without optimization
		double h = rg_->h().minCoeff();

		for (int i = 0; i < momenta_i.rows(); ++i)
		{
			for (int p = 0; p < particleMasses.size(); ++p)
			{
				Matrix3d B;
				B.row(0) = particleAffineMomenta_1.row(p);
				B.row(1) = particleAffineMomenta_2.row(p);
				B.row(2) = particleAffineMomenta_3.row(p);

				Vector3d xp = particlePositions.row(p);
				Vector3d xi = rg_->positions().row(i);

				momenta_i.row(i) += omegas.coeff(p, i) * particleMasses[p]
					* 3.0 / h / h * B * (xi - xp);
			}
		}
		*/

		double h = rg_->h().minCoeff();

		double inertiaTensorRatio = 3.0 / h / h;
		const MatrixX3d& cellAffineMomenta_1 = inertiaTensorRatio
			* omegas.transpose() * particleMasses.asDiagonal()
			* particleAffineMomenta_1;
		const MatrixX3d& cellAffineMomenta_2 = inertiaTensorRatio
			* omegas.transpose() * particleMasses.asDiagonal()
			* particleAffineMomenta_2;
		const MatrixX3d& cellAffineMomenta_3 = inertiaTensorRatio
			* omegas.transpose() * particleMasses.asDiagonal()
			* particleAffineMomenta_3;

		momenta_i.col(0) += (cellAffineMomenta_1.array() * rg_->positions().array())
			.rowwise().sum().matrix();
		momenta_i.col(1) += (cellAffineMomenta_2.array() * rg_->positions().array())
			.rowwise().sum().matrix();
		momenta_i.col(2) += (cellAffineMomenta_3.array() * rg_->positions().array())
			.rowwise().sum().matrix();

		momenta_i.col(0) -= inertiaTensorRatio * omegas.transpose()
			* particleMasses.cwiseProduct(
			(particleAffineMomenta_1.array() * particlePositions.array()).rowwise().sum().matrix());
		momenta_i.col(1) -= inertiaTensorRatio * omegas.transpose()
			* particleMasses.cwiseProduct(
			(particleAffineMomenta_2.array() * particlePositions.array()).rowwise().sum().matrix());
		momenta_i.col(2) -= inertiaTensorRatio * omegas.transpose()
			* particleMasses.cwiseProduct(
			(particleAffineMomenta_3.array() * particlePositions.array()).rowwise().sum().matrix());
	};
	
	if (ps_ != nullptr)
	{
		ADD_AFFINE_MOMENTA(omegas_,
			ps_->masses,
			ps_->affineMomenta_1,
			ps_->affineMomenta_2,
			ps_->affineMomenta_3,
			ps_->positions);
	}

	if (mesh_ != nullptr)
	{
		ADD_AFFINE_MOMENTA(vertexOmegas_,
			mesh_->vertexMasses,
			mesh_->vertexAffineMomenta_1,
			mesh_->vertexAffineMomenta_2,
			mesh_->vertexAffineMomenta_3,
			mesh_->vertexPositions);
		ADD_AFFINE_MOMENTA(elementOmegas_,
			mesh_->elementMasses,
			mesh_->elementAffineMomenta_1,
			mesh_->elementAffineMomenta_2,
			mesh_->elementAffineMomenta_3,
			mesh_->elementPositions);
	}
#endif

	//fout << "momemtum grid after APIC: "
	//	<< momenta_i.colwise().sum() << endl;
	//fout << "angular mometum after APIC: "
	//	<< totalAngularMomentum(rg_->positions(), momenta_i).transpose() << endl;

	rg_->velocities.setZero();
	for (int c = 0; c < rg_->velocities.rows(); ++c)
	{
		if (rg_->masses[c] > 0.0)
		{
			rg_->velocities.row(c) = momenta_i.row(c) / rg_->masses[c];
		}
	}

	if (evaluateVolumesAndDensities)
	{
		if (ps_ != nullptr)
		{
			ps_->densities = omegas_ * rg_->masses / rg_->gridVolume();
			ps_->volumes = ps_->masses.cwiseProduct(ps_->densities.cwiseInverse());
		}
	}
}

void HybridSolver::computeGridForces_(double Dt, MaterialType type)
{
	rg_->forces.setZero();

	if (ps_ != nullptr)
	{
		int numParticles = ps_->elasticDeformationGradients.size();
		int numGrids = rg_->gridNumber();
		candidateElasticDeformationGradients_.resize(numParticles);

		double E0 = ps_->youngsModulus;
		double nu0 = ps_->poissonRatio;

		double lambda0 = E0 * nu0 / (1.0 + nu0) / (1.0 - 2.0 * nu0);
		double mu0 = E0 / 2.0 / (1.0 + nu0);

		double hardeningCoeff = 10.0;

		MatrixX3d EDGModifiers_col_1 = Dt * domegas_1_ * rg_->velocities;
		MatrixX3d EDGModifiers_col_2 = Dt * domegas_2_ * rg_->velocities;
		MatrixX3d EDGModifiers_col_3 = Dt * domegas_3_ * rg_->velocities;

		VectorXd cauchyStress_11_(numParticles), cauchyStress_12_(numParticles), cauchyStress_13_(numParticles);
		VectorXd cauchyStress_21_(numParticles), cauchyStress_22_(numParticles), cauchyStress_23_(numParticles);
		VectorXd cauchyStress_31_(numParticles), cauchyStress_32_(numParticles), cauchyStress_33_(numParticles);

		for (int p = 0; p < numParticles; ++p)
		{
			double lambda, mu;

			if (type == SNOW)
			{
				double J_plastic = ps_->plasticDeformationGradients[p].determinant();

				lambda = lambda0 * exp(hardeningCoeff * (1 - J_plastic));
				mu = mu0 * exp(hardeningCoeff * (1 - J_plastic));
			}
			else if (type == SAND)
			{
				lambda = lambda0;
				mu = mu0;
			}


			const Vector3d& posRef = ps_->positions.row(p).transpose() - rg_->minBound();

			Matrix3d EDG_modifier;
			EDG_modifier.setZero();
			EDG_modifier.col(0) = EDGModifiers_col_1.row(p);
			EDG_modifier.col(1) = EDGModifiers_col_2.row(p);
			EDG_modifier.col(2) = EDGModifiers_col_3.row(p);

			const MatrixX3d& EDG = ps_->elasticDeformationGradients[p];
			const MatrixX3d& PDG = ps_->plasticDeformationGradients[p];

			MatrixX3d virtualEDG = EDG + EDG_modifier * EDG;
			candidateElasticDeformationGradients_[p] = virtualEDG;

			JacobiSVD<Matrix3d> svd(virtualEDG, ComputeFullU | ComputeFullV);

			Matrix3d cauchyStress;

			cauchyStress.setZero();

			if (type == SNOW)
			{
				// polar decomposition
				const Matrix3d& stretchMatrix = svd.matrixV() * svd.singularValues().asDiagonal()
					* svd.matrixV().transpose();
				const Matrix3d& rotationMatrix = svd.matrixU() * svd.matrixV().transpose();

				cauchyStress = ps_->volumes[p] * (2.0 * mu * (virtualEDG - rotationMatrix)
					+ lambda * (virtualEDG.determinant() - 1.0)
					* virtualEDG.determinant() * virtualEDG.transpose().inverse())
					* EDG.transpose();
			}
			else if (type == SAND)
			{
				const Vector3d& Sigma = svd.singularValues();

				Vector3d lnSigma = Vector3d(log(Sigma[0]), log(Sigma[1]), log(Sigma[2]));

				cauchyStress = ps_->volumes[p] * (
					svd.matrixU() * (
						2 * mu * (Sigma.cwiseInverse().cwiseProduct(lnSigma))
						+ lambda * lnSigma.sum() * Sigma.cwiseInverse()
						).asDiagonal()
					* svd.matrixV().transpose()
					) * EDG.transpose();
			}

			cauchyStress_11_[p] = cauchyStress(0, 0);
			cauchyStress_12_[p] = cauchyStress(0, 1);
			cauchyStress_13_[p] = cauchyStress(0, 2);

			cauchyStress_21_[p] = cauchyStress(1, 0);
			cauchyStress_22_[p] = cauchyStress(1, 1);
			cauchyStress_23_[p] = cauchyStress(1, 2);

			cauchyStress_31_[p] = cauchyStress(2, 0);
			cauchyStress_32_[p] = cauchyStress(2, 1);
			cauchyStress_33_[p] = cauchyStress(2, 2);
		}

		

		rg_->forces.col(0) -= domegas_1_.transpose() * cauchyStress_11_;
		rg_->forces.col(0) -= domegas_2_.transpose() * cauchyStress_12_;
		rg_->forces.col(0) -= domegas_3_.transpose() * cauchyStress_13_;

		rg_->forces.col(1) -= domegas_1_.transpose() * cauchyStress_21_;
		rg_->forces.col(1) -= domegas_2_.transpose() * cauchyStress_22_;
		rg_->forces.col(1) -= domegas_3_.transpose() * cauchyStress_23_;

		rg_->forces.col(2) -= domegas_1_.transpose() * cauchyStress_31_;
		rg_->forces.col(2) -= domegas_2_.transpose() * cauchyStress_32_;
		rg_->forces.col(2) -= domegas_3_.transpose() * cauchyStress_33_;

	}

	//clog << "strange negative force: " << rg_->forces.row(482125) << endl;
	if (mesh_ != nullptr)
	{
		// evaluate the in plane forces regarding the mesh
		MatrixX3d vertexInPlaneForces;
		vector<Matrix2d> inPlanePiolaKirhoffStresses;

		mesh_->computeVertexInPlaneForces(vertexInPlaneForces,
			inPlanePiolaKirhoffStresses);
		rg_->forces += vertexOmegas_.transpose() * vertexInPlaneForces;

		// evaluate the normal forces
		//int Nf = mesh_->faces.rows();

		//VectorXd meshStress_11(Nf), meshStress_12(Nf), meshStress_13(Nf),
		//	meshStress_21(Nf), meshStress_22(Nf), meshStress_23(Nf),
		//	meshStress_31(Nf), meshStress_32(Nf), meshStress_33(Nf);

		//for (int f = 0; f < mesh_->faces.rows(); ++f)
		//{
		//	Matrix3d restDirectionMatrix;
		//	restDirectionMatrix.col(0) = mesh_->elementRestDirections_1().row(f);
		//	restDirectionMatrix.col(1) = mesh_->elementRestDirections_2().row(f);
		//	restDirectionMatrix.col(2) = mesh_->elementRestDirections_3().row(f);

		//	Matrix3d directionMatrix;
		//	directionMatrix.col(0) = mesh_->elementDirections_1.row(f);
		//	directionMatrix.col(1) = mesh_->elementDirections_2.row(f);
		//	directionMatrix.col(2) = mesh_->elementDirections_3.row(f);

		//	Matrix3d Q, R;
		//	GramSchmidtOrthonomalization(Q, R, directionMatrix);

		//	const Matrix2d& inPlaneDr = inPlanePiolaKirhoffStresses[f];
		//	double dr11 = inPlaneDr(0, 0);
		//	double dr12 = inPlaneDr(0, 1);
		//	double dr22 = inPlaneDr(1, 1);

		//	double dr13 = mesh_->shearStiffness * R(0, 2);
		//	double dr23 = mesh_->shearStiffness * R(1, 2);

		//	double dr33 = R(2, 2) > 1.0 ? 0.0 : mesh_->stiffness * (1.0 - R(2, 2))
		//		* (1.0 - R(2, 2));

		//	Matrix3d dR;
		//	dR << dr11, dr12, dr13,
		//		0.0, dr22, dr23,
		//		0.0, 0.0, dr33;
		//	Matrix3d K = dR * R.transpose();

		//	Vector3d dF3 = Q * (K.triangularView<StrictlyUpper>().toDenseMatrix()
		//		+ K.triangularView<Upper>().transpose().toDenseMatrix()) * R.inverse().transpose()
		//		* restDirectionMatrix.transpose().col(2);

		//	Matrix3d meshStress = mesh_->elementVolumes[f] * dF3 * directionMatrix.col(2).transpose();

		//	meshStress_11[f] = meshStress(0, 0);
		//	meshStress_12[f] = meshStress(0, 1);
		//	meshStress_13[f] = meshStress(0, 2);

		//	meshStress_21[f] = meshStress(1, 0);
		//	meshStress_22[f] = meshStress(1, 1);
		//	meshStress_23[f] = meshStress(1, 2);

		//	meshStress_31[f] = meshStress(2, 0);
		//	meshStress_32[f] = meshStress(2, 1);
		//	meshStress_33[f] = meshStress(2, 2);
		//}

		//rg_->forces.col(0) -= delementOmegas_1_.transpose() * meshStress_11;
		//rg_->forces.col(0) -= delementOmegas_2_.transpose() * meshStress_12;
		//rg_->forces.col(0) -= delementOmegas_3_.transpose() * meshStress_13;

		//rg_->forces.col(1) -= delementOmegas_1_.transpose() * meshStress_21;
		//rg_->forces.col(1) -= delementOmegas_2_.transpose() * meshStress_22;
		//rg_->forces.col(1) -= delementOmegas_3_.transpose() * meshStress_23;

		//rg_->forces.col(2) -= delementOmegas_1_.transpose() * meshStress_31;
		//rg_->forces.col(2) -= delementOmegas_2_.transpose() * meshStress_32;
		//rg_->forces.col(2) -= delementOmegas_3_.transpose() * meshStress_33;
	}
}

void HybridSolver::gridCollisionHandling_(Eigen::MatrixX3d& gridVelocityChanges)
{
	// consider the mesh and the particle share the same friction
	double friction = 0.2;

	for (int k = 0; k < rg_->resolution()[2]; ++k)
	{
		for (int j = 0; j < rg_->resolution()[1]; ++j)
		{
			for (int i = 0; i < rg_->resolution()[0]; ++i)
			{
				Vector3d gridPos(
					rg_->minBound()[0] + i * rg_->h()[0],
					rg_->minBound()[1] + j * rg_->h()[1],
					rg_->minBound()[2] + k * rg_->h()[2]);

				if (phi_(gridPos) <= 0.0)
				{
					int index = rg_->toIndex(i, j, k);
					const Vector3d& velocity = rg_->velocities.row(index);
					const Vector3d& normal = dphi_(gridPos);

					Vector3d vRef = velocity; // for static object

					double vNormal = vRef.dot(normal);

					if (vNormal < 0.0) // approaching
					{
						Vector3d vTangential = vRef - vNormal * normal;
						if (vTangential.norm() < -friction * vNormal)
						{
							vRef = Vector3d::Zero();
						}
						else
						{
							vRef = vTangential 
								+ friction * vNormal * vTangential / vTangential.norm();
						}
						
						Vector3d vContact = vRef;

						rg_->velocities.row(index) = vContact;
						gridVelocityChanges.row(index) += vContact - velocity;
					}
				}
			}
		}
	}

	if (mesh_ != nullptr)
	{
		for (int s = 0; s < vertexOmegas_.outerSize(); ++s)
			for (SparseMatrix<double>::InnerIterator it(vertexOmegas_, s); it; ++it)
			{
				int vertexID = it.row();
				int gridID = it.col();

				if (mesh_->vertexIsFixed(vertexID))
				{
					int ri = get<0>(rg_->toCoordinate(gridID));
					int rj = get<1>(rg_->toCoordinate(gridID));
					int rk = get<2>(rg_->toCoordinate(gridID));

					for (int i = ri - 1; i <= ri + 1; ++i)
					{
						for (int j = rj - 1; j <= rj + 1; ++j)
						{
							for (int k = rk - 1; k <= rk + 1; ++k)
							{
								if ((i - ri) * (i - ri)
									+ (j - rj) * (j - rj)
									+ (k - rk) * (k - rk)
									<= 0.25)
								{
									int index = rg_->toIndex(i, j, k);
									if (0 <= index && index < rg_->gridNumber())
									{
										rg_->velocities.row(index).setZero();
										gridVelocityChanges.row(index).setZero();
									}
								}
							}
						}
					}
				}
			}
	}
}

void HybridSolver::updateDeformationGradients_(double Dt, MaterialType type)
{
	if (ps_ != nullptr)
	{
		for (int p = 0; p < ps_->elasticDeformationGradients.size(); ++p)
		{
			Matrix3d totalDeformationGradient = candidateElasticDeformationGradients_[p]
				* ps_->plasticDeformationGradients[p];
			JacobiSVD<Matrix3d> svd(candidateElasticDeformationGradients_[p], ComputeFullU | ComputeFullV);

			Matrix3d U = svd.matrixU();
			Vector3d Sigma = svd.singularValues();
			Matrix3d V = svd.matrixV();

			if (type == SNOW)
			{
				Sigma[0] = clamp(Sigma[0], 1.0 - ps_->criticalCompression, 1.0 + ps_->criticalStretch);
				Sigma[1] = clamp(Sigma[1], 1.0 - ps_->criticalCompression, 1.0 + ps_->criticalStretch);
				Sigma[2] = clamp(Sigma[2], 1.0 - ps_->criticalCompression, 1.0 + ps_->criticalStretch);
			}
			else if (type == SAND)
			{
				double E0 = ps_->youngsModulus;
				double nu0 = ps_->poissonRatio;

				double lambda = E0 * nu0 / (1.0 + nu0) / (1.0 - 2.0 * nu0);
				double mu = E0 / 2.0 / (1.0 + nu0);

				// hardening update
				double h0 = 35.0;
				double h1 = 9.0;
				double h2 = 0.2;
				double h3 = 10.0;

				double frictionAngle = (h0 + (h1 * ps_->plasticAmount[p] - h3)
					* exp(-h2 * ps_->plasticAmount[p])) * igl::PI / 180.0;

				double hardeningCoeff = sqrt(2.0 / 3.0) * 2.0 * sin(frictionAngle)
					/ (3.0 - sin(frictionAngle));

				Vector3d lnSigma(log(Sigma[0]), log(Sigma[1]), log(Sigma[2]));
				Vector3d varlnSigma = lnSigma - lnSigma.sum() / 3.0 * Vector3d::Ones();
				double plasticDeformationAmount =
					varlnSigma.norm() + (3.0 * lambda + 2.0 * mu) / 2.0 / mu
					* lnSigma.sum() * hardeningCoeff;

				if (plasticDeformationAmount <= 0.0)
				{
					// do nothing Sigma = Sigma
				}
				else if (varlnSigma.norm() == 0.0 || lnSigma.sum() > 0.0)
				{
					Sigma = Vector3d::Ones();
					ps_->plasticAmount[p] += lnSigma.norm();
				}
				else
				{
					Vector3d Hp = lnSigma
						- plasticDeformationAmount * varlnSigma / varlnSigma.norm();
					Sigma = Vector3d(exp(Hp[0]), exp(Hp[1]), exp(Hp[2]));
					ps_->plasticAmount[p] += plasticDeformationAmount;
				}
			}
			ps_->elasticDeformationGradients[p] = U * Sigma.asDiagonal() * V.transpose();
			ps_->plasticDeformationGradients[p] = V * Sigma.cwiseInverse().asDiagonal()
				* U.transpose() * totalDeformationGradient;


		}
	}


	if (mesh_ != nullptr)
	{
		MatrixX3d candidateGridPositions = rg_->positions() + Dt * rg_->velocities;
		MatrixX3d candidateVertexPositions = mesh_->vertexPositions
			+ Dt * vertexOmegas_ * rg_->velocities;

		const MatrixX3d& elementMultiplier_col_1 = delementOmegas_1_ * candidateGridPositions;
		const MatrixX3d& elementMultiplier_col_2 = delementOmegas_2_ * candidateGridPositions;
		const MatrixX3d& elementMultiplier_col_3 = delementOmegas_3_ * candidateGridPositions;

		for (int f = 0; f < mesh_->elementDirections_1.rows(); ++f)
		{
			Matrix3d restDirectionMatrix;
			restDirectionMatrix.col(0) = mesh_->elementRestDirections_1().row(f);
			restDirectionMatrix.col(1) = mesh_->elementRestDirections_2().row(f);
			restDirectionMatrix.col(2) = mesh_->elementRestDirections_3().row(f);

			Vector3d d1 = candidateVertexPositions.row(mesh_->faces.row(f)[1])
				- candidateVertexPositions.row(mesh_->faces.row(f)[0]);
			Vector3d d2 = candidateVertexPositions.row(mesh_->faces.row(f)[2])
				- candidateVertexPositions.row(mesh_->faces.row(f)[0]);

			Matrix3d elementMul;
			elementMul.col(0) = elementMultiplier_col_1.row(f);
			elementMul.col(1) = elementMultiplier_col_2.row(f);
			elementMul.col(2) = elementMultiplier_col_3.row(f);

			Vector3d d3 = elementMul * (mesh_->elementDirections_3.row(f).transpose());

			Matrix3d virtualElementDirection;
			virtualElementDirection.col(0) = d1;
			virtualElementDirection.col(1) = d2;
			virtualElementDirection.col(2) = d3;

			Matrix3d Q, R;
			GramSchmidtOrthonomalization(Q, R, virtualElementDirection);

			// main return mapping algorithm

			if (R(2, 2) > 1.0)
			{
				R(2, 2) = 1.0;
				R(0, 2) = R(1, 2) = 0.0;
			}
			else
			{
				double normalForce = mesh_->stiffness * (R(2, 2) - 1.0) * (R(2, 2) - 1.0);
				double shearForce = mesh_->shearStiffness * sqrt(
					R(0, 2) * R(0, 2) + R(1, 2) * R(1, 2));


				if (shearForce > mesh_->frictionCoeff * normalForce)
				{
					R(0, 2) *= mesh_->frictionCoeff * normalForce / shearForce;
					R(1, 2) *= mesh_->frictionCoeff * normalForce / shearForce;
				}
			}

			d3 = Q * R.col(2);

			mesh_->elementDirections_1.row(f) = d1;
			mesh_->elementDirections_2.row(f) = d2;
			mesh_->elementDirections_3.row(f) = d3;

		}
	}
}

void HybridSolver::updateGridVelocities_(double Dt, double tolerance,
	Eigen::MatrixX3d& gridVelocityChanges)
{
	for (int c = 0; c < rg_->velocities.rows(); ++c)
	{
		Vector3d force = rg_->forces.row(c);

		if (rg_->masses[c] > 0.0)
		{
			Vector3d vel_add = Dt * force / rg_->masses[c];


			gridVelocityChanges.row(c) = vel_add;

			rg_->velocities.row(c) += vel_add;
		}

		gridVelocityChanges.row(c) += Dt * Vector3d(0.0, 0.0, -9.8);

		rg_->velocities.row(c) += Dt * Vector3d(0.0, 0.0, -9.8);
	}
}

void HybridSolver::updateParticleVelocities_(double alpha, double Dt,
	const Eigen::MatrixX3d& gridVelocityChanges)
{

#if defined(AFFINE_PIC)
	if (ps_ != nullptr)
	{
		ps_->velocities = omegas_ * rg_->velocities;
	}

	if (mesh_ != nullptr)
	{
		mesh_->vertexVelocities = vertexOmegas_ * rg_->velocities;
		for (int f = 0; f < mesh_->elementVelocities.rows(); ++f)
		{
			Vector3d vel_1 = mesh_->vertexVelocities.row(mesh_->faces.row(f)[0]);
			Vector3d vel_2 = mesh_->vertexVelocities.row(mesh_->faces.row(f)[1]);
			Vector3d vel_3 = mesh_->vertexVelocities.row(mesh_->faces.row(f)[2]);

			mesh_->elementVelocities.row(f) = (vel_1 + vel_2 + vel_3) / 3.0;
		}
	}

#else
	if (ps_ != nullptr)
	{
		ps_->velocities = (1.0 - alpha) * omegas_ * rg_->velocities
			+ alpha * (ps_->velocities + Dt * omegas_ * gridVelocityChanges);
	}
	if (mesh_ != nullptr)
	{
		mesh_->elementVelocities = (1.0 - alpha) * elementOmegas_ * rg_->velocities
			+ alpha * (mesh_->elementVelocities + Dt * elementOmegas_ * gridVelocityChanges);
		mesh_->vertexVelocities = (1.0 - alpha) * vertexOmegas_ * rg_->velocities
			+ alpha * (mesh_->vertexVelocities + Dt * vertexOmegas_ * gridVelocityChanges);
	}
#endif
}

void HybridSolver::updateAffineMomenta_(
	Eigen::MatrixX3d& affineMomenta_1,
	Eigen::MatrixX3d& affineMomenta_2,
	Eigen::MatrixX3d& affineMomenta_3,
	const Eigen::SparseMatrix<double>& omegas,
	const Eigen::MatrixX3d& particlePositions,
	const Eigen::MatrixX3d& particleVelocities,
	double dampRatio)
{
	/* the original code without optimization 

	affineMomenta_1.resize(particlePositions.rows(), 3);
	affineMomenta_2.resize(particlePositions.rows(), 3);
	affineMomenta_3.resize(particlePositions.rows(), 3);

	affineMomenta_1.setZero();
	affineMomenta_2.setZero();
	affineMomenta_3.setZero();

	for (int p = 0; p < affineMomenta_1.rows(); ++p)
	{
		Matrix3d B;
		B.setZero();

		for (int i = 0; i < rg_->positions().rows(); ++i)
		{
			Vector3d xp = particlePositions.row(p);
			Vector3d xi = rg_->positions().row(i);
			B += omegas.coeff(p, i) * rg_->velocities.row(i).transpose() * (xi - xp).transpose();
		}

		affineMomenta_1.row(p) = B.row(0);
		affineMomenta_2.row(p) = B.row(1);
		affineMomenta_3.row(p) = B.row(2);
	}
	*/

	MatrixX3d particleVelocities_t = omegas * rg_->velocities;
	affineMomenta_1 = omegas * rg_->velocities.col(0).asDiagonal()
		* rg_->positions()
		- particleVelocities_t.col(0).asDiagonal() * particlePositions;
	affineMomenta_2 = omegas * rg_->velocities.col(1).asDiagonal()
		* rg_->positions()
		- particleVelocities_t.col(1).asDiagonal() * particlePositions;
	affineMomenta_3 = omegas * rg_->velocities.col(2).asDiagonal()
		* rg_->positions()
		- particleVelocities_t.col(2).asDiagonal() * particlePositions;

	for (int i = 0; i < affineMomenta_1.rows(); ++i)
	{
		Matrix3d C;
		C.row(0) = affineMomenta_1.row(i);
		C.row(1) = affineMomenta_2.row(i);
		C.row(2) = affineMomenta_3.row(i);

		Matrix3d skew, symm;
		symm << C(0, 0), 0.5 * (C(0, 1) + C(1, 0)), 0.5 * (C(0, 2) + C(2, 0)),
			0.5 * (C(0, 1) + C(1, 0)), C(1, 1), 0.5 * (C(1, 2) + C(2, 1)),
			0.5 * (C(0, 2) + C(2, 0)), 0.5 * (C(1, 2) + C(2, 1)), C(2, 2);
		skew = C - symm;
		C = skew + (1 - dampRatio) * symm;
		affineMomenta_1.row(i) = C.row(0);
		affineMomenta_2.row(i) = C.row(1);
		affineMomenta_3.row(i) = C.row(2);
	}
}

void HybridSolver::solve(double CFL, double maxt, double alpha)
{
	ofstream fout("debug_log.txt");
	clog << "evaluate the initial weights...";
	if (ps_ != nullptr)
	{
		evaluateInterpolationWeights_(omegas_,
			domegas_1_,
			domegas_2_,
			domegas_3_,
			ps_->positions);
	}
	if (mesh_ != nullptr)
	{
		evaluateInterpolationWeights_(vertexOmegas_,
			dvertexOmegas_1_,
			dvertexOmegas_2_,
			dvertexOmegas_3_,
			mesh_->vertexPositions);
		evaluateInterpolationWeights_(elementOmegas_,
			delementOmegas_1_,
			delementOmegas_2_,
			delementOmegas_3_,
			mesh_->elementPositions);
	}
	clog << "done!\n";

	clog << "first particle to grid...";
	particleToGrid_(1e-20, true);
	clog << "done!\n";

	double t = 0.0;
	clog << "begin to solve...\n";
	while (t <= maxt)
	{
		fout << "mesh, t = " << t << endl;
		double Dt = CFL / max(3e3, rg_->CFL_condition());
		clog << "time: " << t << ", ";
		clog << "delta t: " << Dt << endl;

		clog << "compute grid forces...";
		computeGridForces_(Dt, SAND);

		for (int s = 0; s < vertexOmegas_.outerSize(); ++s)
			for (SparseMatrix<double>::InnerIterator it(vertexOmegas_, s); it; ++it)
			{
				int vertexID = it.row();
				int gridID = it.col();

				if (vertexID == 8)
				{
					fout << "id: " << gridID << ", force: " << rg_->forces.row(gridID) << endl;
				}
			}
		fout << endl;

		clog << "done!\n";

		clog << "update grid velocities...";
		MatrixX3d gridVelocityChanges;
		gridVelocityChanges.resize(rg_->velocities.rows(), 3);

		updateGridVelocities_(Dt, 1e-20, gridVelocityChanges);
		Dt = CFL / max(3e3, rg_->CFL_condition());
		clog << "delta t: " << Dt << endl;
		//fout << "gridvel" << endl << rg_->velocities << endl;
		clog << "done!\n";

		clog << "grid collision handling...";
		gridCollisionHandling_(gridVelocityChanges);
		clog << "done!\n";

		clog << "update deformation gradients...";
		updateDeformationGradients_(Dt, SAND);
		clog << "done!\n";

		clog << "update particle velocities...";
		updateParticleVelocities_(alpha, Dt, gridVelocityChanges);
		fout <<  mesh_->vertexVelocities << endl << endl;
		//fout << "meshvel" << endl << mesh_->vertexVelocities << endl;
		clog << "done!\n";

#if defined(AFFINE_PIC)
		clog << "update affine momenta...";
		if (ps_ != nullptr)
		{
			updateAffineMomenta_(ps_->affineMomenta_1,
				ps_->affineMomenta_2,
				ps_->affineMomenta_3,
				omegas_,
				ps_->positions,
				ps_->velocities,
				0.0);
		}
		if (mesh_ != nullptr)
		{
			//updateAffineMomenta_(mesh_->vertexAffineMomenta_1,
			//	mesh_->vertexAffineMomenta_2,
			//	mesh_->vertexAffineMomenta_3,
			//	vertexOmegas_,
			//	mesh_->vertexPositions,
			//	mesh_->vertexVelocities,
			//	0.0);
			//updateAffineMomenta_(mesh_->elementAffineMomenta_1,
			//	mesh_->elementAffineMomenta_2,
			//	mesh_->elementAffineMomenta_3,
			//	elementOmegas_,
			//	mesh_->elementPositions,
			//	mesh_->elementVelocities,
			//	0.0);
		}
		clog << "done!\n";


#endif
		clog << "update particle positions...";

		if (ps_ != nullptr)
		{
			ps_->positions += Dt * ps_->velocities;
		}
		if (mesh_ != nullptr)
		{
			mesh_->vertexPositions += Dt * mesh_->vertexVelocities;
			//fout << "particle" << endl << ps_->positions - restpPos << endl;
			fout << mesh_->vertexPositions << endl;
			mesh_->updateElementPositions();
		}
		clog << "done!\n";

		clog << "evaluate iteration weights...";
		if (ps_ != nullptr)
		{
			evaluateInterpolationWeights_(omegas_,
				domegas_1_,
				domegas_2_,
				domegas_3_,
				ps_->positions);
		}
		if (mesh_ != nullptr)
		{
			evaluateInterpolationWeights_(vertexOmegas_,
				dvertexOmegas_1_,
				dvertexOmegas_2_,
				dvertexOmegas_3_,
				mesh_->vertexPositions);
			evaluateInterpolationWeights_(elementOmegas_,
				delementOmegas_1_,
				delementOmegas_2_,
				delementOmegas_3_,
				mesh_->elementPositions);
		}
		clog << "done!\n";

		clog << "iterate particle to grid...";
		particleToGrid_(1e-20, false);
		//clog << "strange velocity: " << rg_->velocities.row(482124) << endl;
		clog << "done!\n";

		t += Dt;
	}

}

void HybridSolver::bindViewer(igl::viewer::Viewer * viewer)
{
	viewer_ = viewer;
	if (ps_ != nullptr)
	{
		ps_->bindViewer(viewer);
	}

	if (rg_ != nullptr)
	{
		rg_->bindViewer(viewer);
	}

	if (mesh_ != nullptr)
	{
		mesh_->bindViewer(viewer);
	}
}

void HybridSolver::updateViewer()
{
	using namespace viewer;
	viewer_->data.clear();
	mtx_.lock();
	if (ps_ != nullptr)
	{
		ps_->updateViewer();
	}
	if (mesh_ != nullptr)
	{
		mesh_->updateViewer();
	}
	mtx_.unlock();
}
