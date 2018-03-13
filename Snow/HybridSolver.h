#pragma once
#include <Eigen/Sparse>

class ParticleSystem;
class RegularGrid;

namespace igl
{
	namespace viewer
	{
		class Viewer;
	}
}

class HybridSolver
{
private:
	ParticleSystem * ps_;
	RegularGrid * rg_;

	igl::viewer::Viewer * viewer_;

	Eigen::SparseMatrix<double> omegas_;
	void evaluateInterpolationWeights_();
public:
	HybridSolver(ParticleSystem * ps = nullptr, RegularGrid * rg = nullptr) :
		ps_(ps), rg_(rg), viewer_(nullptr) {}

	void setParticleSystem(ParticleSystem * ps) { ps_ = ps; }
	void setRegularGrid(RegularGrid * rg) { rg_ = rg; }
	void solve(double Dt);

	void bindViewer(igl::viewer::Viewer * viewer);
	void updateViewer();

};