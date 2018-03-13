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
	igl::viewer::Viewer * viewer_;
public:

	Eigen::VectorXd masses;
	Eigen::VectorXd volumes;
	Eigen::VectorXd densities;

	Eigen::MatrixX3d velocities;
	Eigen::MatrixX3d positions;

	std::vector<Eigen::Matrix3d> deformationGradients;
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