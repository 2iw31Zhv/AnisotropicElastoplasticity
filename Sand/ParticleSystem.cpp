#include "ParticleSystem.h"
#include <igl/viewer/Viewer.h>
#include <random>


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
	const Eigen::VectorXd& plasticAmount,
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
	plasticAmount(plasticAmount),
	youngsModulus(youngsModulus),
	poissonRatio(poissonRatio),
	criticalCompression(criticalCompression),
	criticalStretch(critialStretch),
	friction(friction),
	viewer_(nullptr)
{
	int N = masses.size();

	affineMomenta_1.resize(N, 3);
	affineMomenta_2.resize(N, 3);
	affineMomenta_3.resize(N, 3);

	sampleSandColor_();
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
	vector<Matrix3d> affineMomenta;

	for (int i = 0; i < sampleNumber; ++i)
	{
		elasticDeformationGradients.push_back(Matrix3d::Identity());
		plasticDeformationGradients.push_back(Matrix3d::Identity());
		affineMomenta.push_back(Matrix3d::Zero());
	}

	double totalMass = 100.0 * 3.14 * radius * radius * radius;

	VectorXd masses(sampleNumber);
	masses.setOnes();
	masses *= totalMass / sampleNumber;
	VectorXd volumes(sampleNumber);
	volumes.setOnes();
	VectorXd densities(sampleNumber);
	densities.setOnes();

	VectorXd plasticAmount(sampleNumber);
	plasticAmount.setZero();

	return ParticleSystem(velocities,
		positions,
		elasticDeformationGradients,
		plasticDeformationGradients,
		masses,
		volumes,
		densities,
		plasticAmount,
		1.4e5,
		0.2,
		2.5e-2,
		7.5e-3,
		0.2);
}

ParticleSystem ParticleSystem::SandBall(const Eigen::Vector3d & center, double radius, int sampleNumber)
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
	vector<Matrix3d> affineMomenta;

	for (int i = 0; i < sampleNumber; ++i)
	{
		elasticDeformationGradients.push_back(Matrix3d::Identity());
		plasticDeformationGradients.push_back(Matrix3d::Identity());
		affineMomenta.push_back(Matrix3d::Zero());
	}

	double totalMass = 2200 * 3.14 * radius * radius * radius;

	VectorXd masses(sampleNumber);
	masses.setOnes();
	masses *= totalMass / sampleNumber;
	VectorXd volumes(sampleNumber);
	volumes.setOnes();
	VectorXd densities(sampleNumber);
	densities.setOnes();

	VectorXd plasticAmount(sampleNumber);
	plasticAmount.setZero();

	return ParticleSystem(velocities,
		positions,
		elasticDeformationGradients,
		plasticDeformationGradients,
		masses,
		volumes,
		densities,
		plasticAmount,
		3.537e5,
		0.3,
		2.5e-2, // not used 
		7.5e-3, // not used
		0.2); // not used
}

void ParticleSystem::updateViewer()
{
	viewer_->data.add_points(positions, colors_);
}

void ParticleSystem::sampleSandColor_()
{
	Matrix3d colorBasis;
	colorBasis.row(0) = Vector3d(225.0, 169.0, 95.0) / 255.0; //yellow
	colorBasis.row(1) = Vector3d(107.0, 84.0, 30.0) / 255.0; //brown
	colorBasis.row(2) = Vector3d(255.0, 255.0, 255.0) / 255.0; //white

	default_random_engine generator;
	discrete_distribution<int> distribution {0.85, 0.1, 0.05};

	int N = masses.size();
	colors_.resize(N, 3);

	for (int i = 0; i < N; ++i)
	{
		colors_.row(i) = colorBasis.row(distribution(generator));
	}
}
