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

	Eigen::MatrixX3d colors_;
	void sampleSandColor_();
public:

	Eigen::VectorXd masses;
	Eigen::VectorXd volumes;
	Eigen::VectorXd densities;

	Eigen::MatrixX3d velocities;
	Eigen::MatrixX3d positions;

	Eigen::MatrixX3d affineMomenta_1;
	Eigen::MatrixX3d affineMomenta_2;
	Eigen::MatrixX3d affineMomenta_3;
	
	std::vector<Eigen::Matrix3d> elasticDeformationGradients;
	std::vector<Eigen::Matrix3d> plasticDeformationGradients;

	double youngsModulus;
	double poissonRatio;
	double criticalCompression;
	double criticalStretch;
	double friction;

	ParticleSystem(
		const Eigen::MatrixX3d& velocities,
		const Eigen::MatrixX3d& positions,
		const std::vector<Eigen::Matrix3d>& elasticDeformationGradients,
		const std::vector<Eigen::Matrix3d>& plasticDeformationGradients,
		const Eigen::VectorXd& masses,
		const Eigen::VectorXd& volumes,
		const Eigen::VectorXd& densities,
		double youngsModulus,
		double poissonRatio,
		double criticalCompression,
		double criticalStretch,
		double friction);

	static ParticleSystem SnowBall(const Eigen::Vector3d& center,
		double radius,
		int sampleNumber);

	static ParticleSystem SandBall(const Eigen::Vector3d& center,
		double radius,
		int sampleNumber);

	void bindViewer(igl::viewer::Viewer * viewer) { viewer_ = viewer; }
	void updateViewer();
	
};