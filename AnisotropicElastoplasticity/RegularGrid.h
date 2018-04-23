#pragma once
#include <vector>
#include <tuple>
#include <Eigen/Core>
#include <mutex>

namespace igl
{
	namespace viewer
	{
		class Viewer;
	}
}

class RegularGrid
{
private:
	Eigen::VectorXd minBound_;
	Eigen::VectorXd maxBound_;

	Eigen::Vector3i resolution_;
	Eigen::Vector3d h_;
	igl::viewer::Viewer * viewer_;

	Eigen::MatrixX3d points_1_;
	Eigen::MatrixX3d points_2_;
	Eigen::MatrixX3d colors_;

	Eigen::MatrixX3d positions_;

	void initializeRenderingData_();
	
	std::mutex mtx_;
public:
	void recomputeColors();
	Eigen::MatrixX3d velocities;
	Eigen::MatrixX3d forces;
	

	Eigen::VectorXd masses;
	
	RegularGrid(const Eigen::VectorXd& minBound,
		const Eigen::VectorXd& maxBound,
		const Eigen::Vector3i& resolution);

	int toIndex(int i, int j, int k) const;

	std::tuple<int, int, int> toCoordinate(int index) const;

	void bindViewer(igl::viewer::Viewer * viewer) { viewer_ = viewer; }
	void updateViewer();
	int gridNumber() const { return resolution_[0] * resolution_[1] * resolution_[2]; }
	double gridVolume() const { return h_[0] * h_[1] * h_[2]; }
	const Eigen::Vector3d& minBound() const { return minBound_; }
	const Eigen::Vector3d& maxBound() const { return maxBound_; }
	const Eigen::Vector3d& h() const { return h_; }
	const Eigen::Vector3i& resolution() const { return resolution_; }
	const Eigen::MatrixX3d& positions() const { return positions_; }
	const double max_velocity() const;
	const double CFL_condition() const { return max_velocity() / h_.minCoeff(); }
};