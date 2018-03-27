#pragma once
#include <Eigen/Sparse>
#include <functional>
#include <mutex>

class ParticleSystem;
class RegularGrid;

namespace igl
{
	namespace viewer
	{
		class Viewer;
	}
}

using LevelSet = std::function<double(const Eigen::Vector3d&)>;
using DLevelSet = std::function<Eigen::Vector3d(const Eigen::Vector3d&)>;

class HybridSolver
{
private:
	ParticleSystem * ps_;
	RegularGrid * rg_;

	igl::viewer::Viewer * viewer_;

	Eigen::SparseMatrix<double> omegas_;
	Eigen::SparseMatrix<double> domegas_1_;
	Eigen::SparseMatrix<double> domegas_2_;
	Eigen::SparseMatrix<double> domegas_3_;

	std::vector<Eigen::Matrix3d> candidateElasticDeformationGradients_;
	
	LevelSet phi_;
	DLevelSet dphi_;


	void evaluateInterpolationWeights_();
	void particleToGrid_();
	void computeGridForces_(double Dt);
	void gridCollisionHandling_();
	void updateDeformationGradients_();
	void updateParticleVelocities_(double alpha, double Dt);

	void updateAffineMomenta_();
public:

	std::mutex mtx_;
	HybridSolver(ParticleSystem * ps = nullptr, RegularGrid * rg = nullptr) :
		ps_(ps), rg_(rg), viewer_(nullptr) {}

	void setParticleSystem(ParticleSystem * ps) { ps_ = ps; }
	void setRegularGrid(RegularGrid * rg) { rg_ = rg; }
	void setLevelSet(const LevelSet& phi, const DLevelSet& dphi) { phi_ = phi; dphi_ = dphi; }
	void solve(double Dt, double maxt, double alpha);

	void bindViewer(igl::viewer::Viewer * viewer);
	void updateViewer();

};