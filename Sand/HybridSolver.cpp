#include "HybridSolver.h"
#include "ParticleSystem.h"
#include "RegularGrid.h"
#include <vector>
#include <igl/viewer/Viewer.h>
#include "interpolation.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace igl;

void HybridSolver::evaluateInterpolationWeights_()
{
	vector<Triplet<double>> omegaTriplets;
	vector<Triplet<double>> domegaTriplets_1;
	vector<Triplet<double>> domegaTriplets_2;
	vector<Triplet<double>> domegaTriplets_3;

	for (int p = 0; p < ps_->positions.rows(); ++p)
	{
		const Vector3d& position = ps_->positions.row(p);
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

						if (dyadic_i > 0 && dyadic_j > 0 && dyadic_k > 0)
						{
							omegaTriplets.push_back(Triplet<double>(p,
								rg_->toIndex(i, j, k),
								dyadic_i * dyadic_j * dyadic_k));
							domegaTriplets_1.push_back(Triplet<double>(p,
								rg_->toIndex(i, j, k),
								ddyadic_i * dyadic_j * dyadic_k));
							domegaTriplets_2.push_back(Triplet<double>(p,
								rg_->toIndex(i, j, k),
								dyadic_i * ddyadic_j * dyadic_k));
							domegaTriplets_3.push_back(Triplet<double>(p,
								rg_->toIndex(i, j, k),
								dyadic_i * dyadic_j * ddyadic_k));
						}
					}
				}
			}
		}
	}

	omegas_.resize(ps_->masses.size(), rg_->masses.size());
	domegas_1_.resize(ps_->masses.size(), rg_->masses.size());
	domegas_2_.resize(ps_->masses.size(), rg_->masses.size());
	domegas_3_.resize(ps_->masses.size(), rg_->masses.size());

	omegas_.setZero();
	domegas_1_.setZero();
	domegas_2_.setZero();
	domegas_3_.setZero();

	omegas_.setFromTriplets(omegaTriplets.begin(), omegaTriplets.end());
	domegas_1_.setFromTriplets(domegaTriplets_1.begin(), domegaTriplets_1.end());
	domegas_2_.setFromTriplets(domegaTriplets_2.begin(), domegaTriplets_2.end());
	domegas_3_.setFromTriplets(domegaTriplets_3.begin(), domegaTriplets_3.end());

}

void HybridSolver::particleToGrid_()
{
	rg_->masses = omegas_.transpose() * ps_->masses;
	rg_->recomputeColors();

	const MatrixX3d& momenta_p = ps_->masses.asDiagonal() * ps_->velocities;
	const MatrixX3d& momenta_i = omegas_.transpose() * momenta_p;
	rg_->velocities = rg_->masses.asDiagonal().inverse() * momenta_i;
}

// need optimization
void HybridSolver::computeGridForces_(double Dt)
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
		double J_plastic = ps_->plasticDeformationGradients[p].determinant();

		double lambda = lambda0 * exp(hardeningCoeff * (1 - J_plastic));
		double mu = mu0 * exp(hardeningCoeff * (1 - J_plastic));

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
		
		// polar decomposition
		JacobiSVD<Matrix3d> svd(virtualEDG, ComputeFullU | ComputeFullV);
		Matrix3d stretchMatrix = svd.matrixV() * svd.singularValues().asDiagonal()
			* svd.matrixV().transpose();
		Matrix3d rotationMatrix = svd.matrixU() * svd.matrixV().transpose();

		Matrix3d cauchyStress = ps_->volumes[p] * (2.0 * mu * (virtualEDG - rotationMatrix)
			+ lambda * (virtualEDG.determinant() - 1.0)
			* virtualEDG.determinant() * virtualEDG.transpose().inverse())
			* EDG.transpose();

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

	for (int c = 0; c < numGrids; ++c)
	{
		rg_->forces.row(c) += rg_->masses[c] * Vector3d(0.0, 0.0, -9.8);
	}
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

}

void HybridSolver::updateDeformationGradients_()
{
	for (int p = 0; p < ps_->elasticDeformationGradients.size(); ++p)
	{
		Matrix3d totalDeformationGradient = candidateElasticDeformationGradients_[p]
			* ps_->plasticDeformationGradients[p];
		JacobiSVD<Matrix3d> svd(candidateElasticDeformationGradients_[p], ComputeFullU | ComputeFullV);

		Matrix3d U = svd.matrixU();
		Vector3d Sigma = svd.singularValues();
		Matrix3d V = svd.matrixV();
		
		Sigma[0] = clamp(Sigma[0], 1.0 - ps_->criticalCompression, 1.0 + ps_->criticalStretch);
		Sigma[1] = clamp(Sigma[1], 1.0 - ps_->criticalCompression, 1.0 + ps_->criticalStretch);
		Sigma[2] = clamp(Sigma[2], 1.0 - ps_->criticalCompression, 1.0 + ps_->criticalStretch);

		ps_->elasticDeformationGradients[p] = U * Sigma.asDiagonal() * V.transpose();
		ps_->plasticDeformationGradients[p] = V * Sigma.cwiseInverse().asDiagonal()
			* U.transpose() * totalDeformationGradient;

	}
}

void HybridSolver::updateParticleVelocities_(double alpha, double Dt)
{
	ps_->velocities = (1.0 - alpha) * omegas_ * rg_->velocities
		+ alpha * (ps_->velocities + Dt * omegas_ * rg_->masses.asDiagonal() * rg_->forces);
}

void HybridSolver::solve(double Dt, double maxt, double alpha)
{
	clog << "evaluate the initial weights...";
	evaluateInterpolationWeights_();
	clog << "done!\n";
	clog << "first particle to grid...";
	particleToGrid_();
	clog << "done!\n";

	clog << "evaluate particle densities and volumes...";
	ps_->densities = omegas_ * rg_->masses / rg_->gridVolume();
	ps_->volumes = ps_->masses.cwiseProduct(ps_->densities.cwiseInverse());
	clog << "done!\n";

	double t = 0.0;

	clog << "begin to solve...\n";
	while (t <= maxt)
	{
		clog << "time: " << t << endl;
		clog << "compute grid forces...";
		computeGridForces_(Dt);
		clog << "done!\n";
		clog << "update grid velocities...";
		rg_->velocities += Dt * rg_->masses.asDiagonal() * rg_->forces;
		clog << "done!\n";
		clog << "grid collision handling...";
		gridCollisionHandling_();
		clog << "done!\n";
		clog << "update deformation gradients...";
		updateDeformationGradients_();
		clog << "done!\n";
		clog << "update particle velocities...";
		updateParticleVelocities_(alpha, Dt);
		clog << "done!\n";
		clog << "update particle positions...";
		ps_->positions += Dt * ps_->velocities;
		clog << "done!\n";
		//updateViewer();

		clog << "evaluate iteration weights...";
		evaluateInterpolationWeights_();
		clog << "done!\n";
		clog << "iterate particle to grid...";
		particleToGrid_();
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
}

void HybridSolver::updateViewer()
{
	using namespace viewer;
	viewer_->data.clear();
	mtx_.lock();
	ps_->updateViewer();
	//rg_->updateViewer();
	MatrixX3d e0, e1, c;
	e0.resize(4, 3);
	e1.resize(4, 3);
	c.resize(4, 3);

	c.setZero();
	c.col(1).setOnes();

	e0.row(0) = Vector3d(-10.0, -10.0, -1.8);
	e0.row(1) = Vector3d(-10.0, 10.0, -1.8);
	e0.row(2) = Vector3d(10.0, 10.0, -1.8);
	e0.row(3) = Vector3d(10.0, -10.0, -1.8);

	e1.row(0) = e0.row(1);
	e1.row(1) = e0.row(2);
	e1.row(2) = e0.row(3);
	e1.row(3) = e0.row(0);

	viewer_->data.add_edges(e0, e1, c);
	mtx_.unlock();
}