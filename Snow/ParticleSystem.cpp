#include "ParticleSystem.h"
#include <igl/viewer/Viewer.h>


using namespace std;
using namespace Eigen;
using namespace igl;

ParticleSystem::ParticleSystem(const Eigen::MatrixX3d& velocities,
	const Eigen::MatrixX3d& positions,
	const std::vector<Eigen::Matrix3d>& deformationGradients,
	const Eigen::VectorXd& masses,
	const Eigen::VectorXd& volumes,
	const Eigen::VectorXd& densities) :
	velocities(velocities),
	positions(positions),
	deformationGradients(deformationGradients),
	masses(masses),
	volumes(volumes),
	densities(densities),
	viewer_(nullptr)
{

}

ParticleSystem ParticleSystem::Ball(const Eigen::Vector3d & center, 
	double radius, 
	int sampleNumber)
{
	MatrixX3d positions;
	positions.resize(sampleNumber, 3);

	for (int i = 0; i < sampleNumber; ++i)
	{
		Vector3d tempPos;
		while ((tempPos = Vector3d::Random()).norm() > 1.0);

		positions.row(i) = center + radius * tempPos;
	}

	MatrixX3d velocities;
	velocities.resize(sampleNumber, 3);
	velocities.setZero();

	vector<Matrix3d> deformationGradients;

	for (int i = 0; i < sampleNumber; ++i)
	{
		deformationGradients.push_back(Matrix3d::Identity());
	}

	VectorXd masses(sampleNumber);
	masses.setOnes();
	VectorXd volumes(sampleNumber);
	volumes.setOnes();
	VectorXd densities(sampleNumber);
	densities.setOnes();

	return ParticleSystem(velocities,
		positions,
		deformationGradients,
		masses,
		volumes,
		densities);
}

void ParticleSystem::updateViewer()
{
	int numParticles = positions.rows();
	using namespace viewer;
	MatrixX3d colors;
	colors.resize(numParticles, 3);
	colors.setOnes();

	viewer_->data.add_points(positions, colors);
}
