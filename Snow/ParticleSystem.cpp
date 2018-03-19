#include "ParticleSystem.h"
#include <igl/viewer/Viewer.h>


using namespace std;
using namespace Eigen;
using namespace igl;

ParticleSystem::ParticleSystem(const Eigen::MatrixX3d& velocities,
	const Eigen::MatrixX3d& positions,
	const std::vector<Eigen::Matrix3d>& elasticDeformationGradients,
	const std::vector<Eigen::Matrix3d>& plasticDeformationGradients,
	const Eigen::VectorXd& masses,
	const Eigen::VectorXd& volumes,
	const Eigen::VectorXd& densities,
	double youngsModulus,
	double poissonRatio,
	double criticalCompression,
	double critialStretch,
	double friction) :
	velocities(velocities),
	positions(positions),
	elasticDeformationGradients(elasticDeformationGradients),
	plasticDeformationGradients(plasticDeformationGradients),
	masses(masses),
	volumes(volumes),
	densities(densities),
	youngsModulus(youngsModulus),
	poissonRatio(poissonRatio),
	criticalCompression(criticalCompression),
	criticalStretch(critialStretch),
	friction(friction),
	viewer_(nullptr)
{

}

ParticleSystem ParticleSystem::SnowBall(const Eigen::Vector3d & center, 
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

	vector<Matrix3d> elasticDeformationGradients;
	vector<Matrix3d> plasticDeformationGradients;

	for (int i = 0; i < sampleNumber; ++i)
	{
		elasticDeformationGradients.push_back(Matrix3d::Identity());
		plasticDeformationGradients.push_back(Matrix3d::Identity());
	}

	VectorXd masses(sampleNumber);
	masses.setOnes();
	VectorXd volumes(sampleNumber);
	volumes.setOnes();
	VectorXd densities(sampleNumber);
	densities.setOnes();

	return ParticleSystem(velocities,
		positions,
		elasticDeformationGradients,
		plasticDeformationGradients,
		masses,
		volumes,
		densities,
		1.4e5,
		0.2,
		2.5e-2,
		7.5e-3,
		0.2);
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
