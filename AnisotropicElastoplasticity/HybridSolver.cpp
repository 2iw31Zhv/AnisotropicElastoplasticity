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

		for (int i = i_fl - 1; i <= i_fl + 2; ++i)
		{
			for (int j = j_fl - 1; j <= j_fl + 2; ++j)
			{
				for (int k = k_fl - 1; k <= k_fl + 2; ++k)
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

}

void HybridSolver::particleToGrid_(
	Eigen::VectorXd& gridMasses,
	Eigen::MatrixX3d& gridVelocities,
	const Eigen::VectorXd& particleMasses,
	const Eigen::MatrixX3d& particleVelocities,
	const Eigen::MatrixX3d& particlePositions,
	const Eigen::MatrixX3d& particleAffineMomenta_1,
	const Eigen::MatrixX3d& particleAffineMomenta_2,
	const Eigen::MatrixX3d& particleAffineMomenta_3,
	const Eigen::SparseMatrix<double>& omegas)
{
	gridMasses += omegas.transpose() * particleMasses;

	const MatrixX3d& momenta_p = particleMasses.asDiagonal() * particleVelocities;
	MatrixX3d momenta_i = omegas.transpose() * momenta_p;

#if defined(AFFINE_PIC)
	double inertiaTensorRatio = 3.0 / pow(rg_->h().prod(), 2.0 / 3.0);
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

#endif
	
	gridVelocities += gridMasses.asDiagonal().inverse() * momenta_i;
}

void HybridSolver::computeGridForces_(double Dt, MaterialType type)
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

	rg_->forces.setZero();

	rg_->forces.col(0) -= domegas_1_.transpose() * cauchyStress_11_;
	rg_->forces.col(0) -= domegas_2_.transpose() * cauchyStress_12_;
	rg_->forces.col(0) -= domegas_3_.transpose() * cauchyStress_13_;

	rg_->forces.col(1) -= domegas_1_.transpose() * cauchyStress_21_;
	rg_->forces.col(1) -= domegas_2_.transpose() * cauchyStress_22_;
	rg_->forces.col(1) -= domegas_3_.transpose() * cauchyStress_23_;

	rg_->forces.col(2) -= domegas_1_.transpose() * cauchyStress_31_;
	rg_->forces.col(2) -= domegas_2_.transpose() * cauchyStress_32_;
	rg_->forces.col(2) -= domegas_3_.transpose() * cauchyStress_33_;


	// evaluate the in plane forces regarding the mesh
	MatrixX3d vertexInPlaneForces;
	vector<Matrix2d> inPlanePiolaKirhoffStresses;

	mesh_->computeVertexInPlaneForces(vertexInPlaneForces,
		inPlanePiolaKirhoffStresses);

	//ofstream fout("inplaneforce.txt", ios::app);
	//fout << endl << vertexInPlaneForces << endl;

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

	//	Matrix3d directionMatrix = mesh_->elasticDeformationGradient[f] * restDirectionMatrix;

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

void HybridSolver::gridCollisionHandling_()
{
	double friction = ps_->friction;

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
					}
				}
			}
		}
	}

	// handle fixed cloth vertex
	for (int s = 0; s < vertexOmegas_.outerSize(); ++s)
		for (SparseMatrix<double>::InnerIterator it(vertexOmegas_, s); it; ++it)
		{
			int vertexID = it.row();
			int gridID = it.col();

			if (vertexID == 0 || vertexID == 24 ||
				vertexID == 600 || vertexID == 624)
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
								<= 1.0)
							{
								int index = rg_->toIndex(i, j, k);
								if (0 <= index && index < rg_->gridNumber())
								{
									rg_->velocities.row(index).setZero();
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

	MatrixX3d candidateVertexPositions = mesh_->vertexPositions
		+ Dt * mesh_->vertexVelocities;

	MatrixX3d elementEDGModifiers_col_1 = Dt * delementOmegas_1_ * rg_->velocities;
	MatrixX3d elementEDGModifiers_col_2 = Dt * delementOmegas_2_ * rg_->velocities;
	MatrixX3d elementEDGModifiers_col_3 = Dt * delementOmegas_3_ * rg_->velocities;

	for (int f = 0; f < mesh_->elasticDeformationGradient.size(); ++f)
	{
		Matrix3d restDirectionMatrix;
		restDirectionMatrix.col(0) = mesh_->elementRestDirections_1().row(f);
		restDirectionMatrix.col(1) = mesh_->elementRestDirections_2().row(f);
		restDirectionMatrix.col(2) = mesh_->elementRestDirections_3().row(f);

		Vector3d d1 = candidateVertexPositions.row(mesh_->faces.row(f)[1])
			- candidateVertexPositions.row(mesh_->faces.row(f)[0]);
		Vector3d d2 = candidateVertexPositions.row(mesh_->faces.row(f)[2])
			- candidateVertexPositions.row(mesh_->faces.row(f)[0]);

		Matrix3d elementMod;
		elementMod.col(0) = elementEDGModifiers_col_1.row(f);
		elementMod.col(1) = elementEDGModifiers_col_2.row(f);
		elementMod.col(2) = elementEDGModifiers_col_3.row(f);

		Vector3d origin_d3 = mesh_->elasticDeformationGradient[f]
			* restDirectionMatrix.col(2);
		Vector3d d3 = origin_d3 + elementMod * origin_d3;

		Matrix3d virtualElementDirection;
		virtualElementDirection.col(0) = d1;
		virtualElementDirection.col(1) = d2;
		virtualElementDirection.col(2) = d3;

		Matrix3d virtualEDG = virtualElementDirection * restDirectionMatrix.inverse();

		Matrix3d totalEDG = virtualEDG * mesh_->plasticDeformationGradient[f];

		Matrix3d Q, R;
		GramSchmidtOrthonomalization(Q, R, virtualElementDirection);

		// main return mapping algorithm

		//if (R(2, 2) > 1.0)
		//{
		//	R(2, 2) = 1.0;
		//	R(0, 2) = R(1, 2) = 0.0;
		//}
		//else
		//{
		//	Vector3d traction = (
		//		mesh_->shearStiffness * (R(0, 2) * Q.col(0)
		//			+ R(1, 2) * Q.col(1))
		//		+ mesh_->stiffness * (R(2, 2) - 1.0) * (R(2, 2) - 1.0) * Q.col(2)
		//		) / R(0, 0) / R(1, 1);
		//	Vector3d n = Q.col(2);
		//	double normalForce = mesh_->stiffness * (R(2, 2) - 1.0) * (R(2, 2) - 1.0);
		//	double shearForce = mesh_->shearStiffness * sqrt(
		//		R(0, 2) * R(0, 2) + R(1, 2) * R(1, 2));


		//	if (shearForce > mesh_->frictionCoeff * normalForce)
		//	{
		//		R(0, 2) *= mesh_->frictionCoeff * normalForce / shearForce;
		//		R(1, 2) *= mesh_->frictionCoeff * normalForce / shearForce;
		//	}
		//}

		mesh_->elasticDeformationGradient[f] = Q * R * restDirectionMatrix.inverse();
		mesh_->plasticDeformationGradient[f] = mesh_->elasticDeformationGradient[f].inverse()
			* totalEDG;

	}
}

void HybridSolver::updateParticleVelocities_(double alpha, double Dt)
{
#if defined(AFFINE_PIC)
	ps_->velocities = omegas_ * rg_->velocities;
	mesh_->elementVelocities = elementOmegas_ * rg_->velocities;
	mesh_->vertexVelocities = vertexOmegas_ * rg_->velocities;

#else
	ps_->velocities = (1.0 - alpha) * omegas_ * rg_->velocities
		+ alpha * (ps_->velocities + Dt * omegas_ * rg_->masses.asDiagonal() * rg_->forces);
	// TODO
#endif
}

void HybridSolver::updateAffineMomenta_(
Eigen::MatrixX3d& affineMomenta_1,
Eigen::MatrixX3d& affineMomenta_2,
Eigen::MatrixX3d& affineMomenta_3,
const Eigen::SparseMatrix<double>& omegas,
const Eigen::MatrixX3d& particlePositions,
const Eigen::MatrixX3d& particleVelocities)
{
	affineMomenta_1 = omegas * rg_->velocities.col(0).asDiagonal()
		* rg_->positions()
		-particleVelocities.col(0).asDiagonal() * particlePositions;
	affineMomenta_2 = omegas * rg_->velocities.col(1).asDiagonal()
		* rg_->positions()
		- particleVelocities.col(1).asDiagonal() * particlePositions;
	affineMomenta_3 = omegas * rg_->velocities.col(2).asDiagonal()
		* rg_->positions()
		- particleVelocities.col(2).asDiagonal() * particlePositions;
}

void HybridSolver::solve(double CFL, double maxt, double alpha)
{
	clog << "evaluate the initial weights...";
	evaluateInterpolationWeights_(omegas_,
		domegas_1_,
		domegas_2_,
		domegas_3_,
		ps_->positions);
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
	clog << "done!\n";

	clog << "first particle to grid...";
	rg_->masses.setOnes();
	rg_->masses *= 1e-20;

	rg_->velocities.setZero();

	particleToGrid_(rg_->masses,
		rg_->velocities,
		ps_->masses,
		ps_->velocities,
		ps_->positions,
		ps_->affineMomenta_1,
		ps_->affineMomenta_2,
		ps_->affineMomenta_3,
		omegas_);

	// DANGEROUS!
	clog << "evaluate particle densities and volumes...";
	ps_->densities = omegas_ * rg_->masses / rg_->gridVolume();
	ps_->volumes = ps_->masses.cwiseProduct(ps_->densities.cwiseInverse());
	clog << "done!\n";

	particleToGrid_(rg_->masses,
		rg_->velocities,
		mesh_->vertexMasses,
		mesh_->vertexVelocities,
		mesh_->vertexPositions,
		mesh_->vertexAffineMomenta_1,
		mesh_->vertexAffineMomenta_2,
		mesh_->vertexAffineMomenta_3,
		vertexOmegas_);
	particleToGrid_(rg_->masses,
		rg_->velocities,
		mesh_->elementMasses,
		mesh_->elementVelocities,
		mesh_->elementPositions,
		mesh_->elementAffineMomenta_1,
		mesh_->elementAffineMomenta_2,
		mesh_->elementAffineMomenta_3,
		elementOmegas_);
	clog << "done!\n";

	double t = 0.0;
	clog << "begin to solve...\n";
	while (t <= maxt)
	{
		double Dt = 0.3 / max(3e3, rg_->CFL_condition());

		clog << "time: " << t << ", ";
		clog << "delta t: " << Dt << endl;

		clog << "compute grid forces...";
		computeGridForces_(Dt, SAND);
		clog << "done!\n";

		
		clog << "update grid velocities...";
		rg_->velocities += Dt * rg_->masses.cwiseInverse().asDiagonal() * rg_->forces;
		
		for (int c = 0; c < rg_->velocities.rows(); ++c)
		{
			rg_->velocities.row(c) += Dt * Vector3d(0.0, 0.0, -9.8);
		}

		clog << "done!\n";

		clog << "grid collision handling...";
		gridCollisionHandling_();
		clog << "done!\n";

		clog << "update deformation gradients...";
		updateDeformationGradients_(Dt, SAND);
		clog << "done!\n";

		clog << "update particle velocities...";
		updateParticleVelocities_(alpha, Dt);
		clog << "done!\n";

#if defined(AFFINE_PIC)
		clog << "update affine momenta...";
		updateAffineMomenta_(ps_->affineMomenta_1,
			ps_->affineMomenta_2,
			ps_->affineMomenta_3,
			omegas_,
			ps_->positions,
			ps_->velocities);
		updateAffineMomenta_(mesh_->vertexAffineMomenta_1,
			mesh_->vertexAffineMomenta_2,
			mesh_->vertexAffineMomenta_3,
			vertexOmegas_,
			mesh_->vertexPositions,
			mesh_->vertexVelocities);
		updateAffineMomenta_(mesh_->elementAffineMomenta_1,
			mesh_->elementAffineMomenta_2,
			mesh_->elementAffineMomenta_3,
			elementOmegas_,
			mesh_->elementPositions,
			mesh_->elementVelocities);
		clog << "done!\n";


#endif
		clog << "update particle positions...";

		ps_->positions += Dt * ps_->velocities;
		mesh_->vertexPositions += Dt * mesh_->vertexVelocities;

		mesh_->updateElementPositions();

		clog << "done!\n";

		clog << "evaluate iteration weights...";
		evaluateInterpolationWeights_(omegas_,
			domegas_1_,
			domegas_2_,
			domegas_3_,
			ps_->positions);
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
		clog << "done!\n";

		clog << "iterate particle to grid...";

		rg_->masses.setOnes();
		rg_->masses *= 1e-20;
		rg_->velocities.setZero();

		particleToGrid_(rg_->masses,
			rg_->velocities,
			ps_->masses,
			ps_->velocities,
			ps_->positions,
			ps_->affineMomenta_1,
			ps_->affineMomenta_2,
			ps_->affineMomenta_3,
			omegas_);
		particleToGrid_(rg_->masses,
			rg_->velocities,
			mesh_->vertexMasses,
			mesh_->vertexVelocities,
			mesh_->vertexPositions,
			mesh_->vertexAffineMomenta_1,
			mesh_->vertexAffineMomenta_2,
			mesh_->vertexAffineMomenta_3,
			vertexOmegas_);
		particleToGrid_(rg_->masses,
			rg_->velocities,
			mesh_->elementMasses,
			mesh_->elementVelocities,
			mesh_->elementPositions,
			mesh_->elementAffineMomenta_1,
			mesh_->elementAffineMomenta_2,
			mesh_->elementAffineMomenta_3,
			elementOmegas_);
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
	ps_->updateViewer();
	mesh_->updateViewer();
	mtx_.unlock();
}
