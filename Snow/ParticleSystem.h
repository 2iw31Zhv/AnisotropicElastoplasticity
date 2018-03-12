#pragma once
#include <vector>
#include <Eigen/Core>

namespace igl
{
	namespace viewer
	{
		class Viewer;
	}
}

class ParticleSystem
{
private:
	
	Eigen::VectorXd masses_;
	Eigen::VectorXd volumes_;
	Eigen::VectorXd densities_;

	Eigen::MatrixX3d velocities_;
	Eigen::MatrixX3d positions_;

	std::vector<Eigen::Matrix3d> deformationGradients_;

	igl::viewer::Viewer * viewer_;
public:
	ParticleSystem(
		const Eigen::MatrixX3d& velocities,
		const Eigen::MatrixX3d& positions,
		const std::vector<Eigen::Matrix3d>& deformationGradients,
		const Eigen::VectorXd& masses,
		const Eigen::VectorXd& volumes,
		const Eigen::VectorXd& densities);

	static ParticleSystem Ball(const Eigen::Vector3d& center,
		double radius,
		int sampleNumber);

	void bindViewer(igl::viewer::Viewer * viewer) { viewer_ = viewer; }
	void updateViewer();
};