#pragma once
#include <Eigen/Sparse>
#include <functional>
#include <mutex>

class ParticleSystem;
class RegularGrid;
class LagrangianMesh;

namespace igl
{
	namespace viewer
	{
		class Viewer;
	}
}

using LevelSet = std::function<double(const Eigen::Vector3d&)>;
using DLevelSet = std::function<Eigen::Vector3d(const Eigen::Vector3d&)>;

enum MaterialType
{
	SNOW = 0,
	SAND
};

class HybridSolver
{
private:
	ParticleSystem * ps_;
	RegularGrid * rg_;
	LagrangianMesh * mesh_;

	igl::viewer::Viewer * viewer_;

	Eigen::SparseMatrix<double> omegas_;
	Eigen::SparseMatrix<double> domegas_1_;
	Eigen::SparseMatrix<double> domegas_2_;
	Eigen::SparseMatrix<double> domegas_3_;

	Eigen::SparseMatrix<double> vertexOmegas_;
	Eigen::SparseMatrix<double> dvertexOmegas_1_;
	Eigen::SparseMatrix<double> dvertexOmegas_2_;
	Eigen::SparseMatrix<double> dvertexOmegas_3_;

	Eigen::SparseMatrix<double> elementOmegas_;
	Eigen::SparseMatrix<double> delementOmegas_1_;
	Eigen::SparseMatrix<double> delementOmegas_2_;
	Eigen::SparseMatrix<double> delementOmegas_3_;


	std::vector<Eigen::Matrix3d> candidateElasticDeformationGradients_;
	
	LevelSet phi_;
	DLevelSet dphi_;


	void evaluateInterpolationWeights_(Eigen::SparseMatrix<double>& omegas,
		Eigen::SparseMatrix<double>& domegas_1,
		Eigen::SparseMatrix<double>& domegas_2,
		Eigen::SparseMatrix<double>& domegas_3,
		const Eigen::MatrixX3d& particlePositions);
	void particleToGrid_(double tolerance, bool evaluateVolumesAndDensities);

	void computeGridForces_(double Dt, MaterialType type);
	void gridCollisionHandling_(Eigen::MatrixX3d& gridVelocityChanges);
	void updateDeformationGradients_(double Dt, MaterialType type);
	void updateGridVelocities_(double Dt, double tolerance,
		Eigen::MatrixX3d& gridVelocityChanges);
	void updateParticleVelocities_(double alpha, double Dt,
		const Eigen::MatrixX3d& gridVelocityChanges);

	void updateAffineMomenta_(Eigen::MatrixX3d& affineMomenta_1,
		Eigen::MatrixX3d& affineMomenta_2,
		Eigen::MatrixX3d& affineMomenta_3,
		const Eigen::SparseMatrix<double>& omegas,
		const Eigen::MatrixX3d& particlePositions,
		const Eigen::MatrixX3d& particleVelocities,
		double dampRatio);
public:

	std::mutex mtx_;
	HybridSolver(ParticleSystem * ps = nullptr, RegularGrid * rg = nullptr) :
		ps_(ps), rg_(rg), viewer_(nullptr) {}

	void setParticleSystem(ParticleSystem * ps) { ps_ = ps; }
	void setRegularGrid(RegularGrid * rg) { rg_ = rg; }
	void setLagrangianMesh(LagrangianMesh * mesh) { mesh_ = mesh; }
	void setLevelSet(const LevelSet& phi, const DLevelSet& dphi) { phi_ = phi; dphi_ = dphi; }
	void solve(double CFL, double maxt, double alpha);

	void bindViewer(igl::viewer::Viewer * viewer);
	void updateViewer();

};