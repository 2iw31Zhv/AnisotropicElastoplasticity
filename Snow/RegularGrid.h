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

struct GridValue
{
	Eigen::Vector3d velocity;
	Eigen::Vector3d force;
	double mass;
};

class RegularGrid
{
private:
	Eigen::VectorXd minBound_;
	Eigen::VectorXd maxBound_;

	Eigen::Vector3i resolution_;
	Eigen::Vector3d h_;
	
	std::vector<GridValue> data_;

	igl::viewer::Viewer * viewer_;

	int toIndex_(int i, int j, int k) const;

	Eigen::MatrixX3d points_1_;
	Eigen::MatrixX3d points_2_;
	Eigen::MatrixX3d colors_;


	void initializeRenderingData_();
public:
	RegularGrid(const Eigen::VectorXd& minBound,
		const Eigen::VectorXd& maxBound,
		const Eigen::Vector3i& resolution);

	const GridValue& at(int i, int j, int k) const;
	void set(int i, int j, int k, const GridValue& gridvalue);

	void bindViewer(igl::viewer::Viewer * viewer) { viewer_ = viewer; }
	void updateViewer();
	int gridNumber() const { return resolution_[0] * resolution_[1] * resolution_[2]; }
};