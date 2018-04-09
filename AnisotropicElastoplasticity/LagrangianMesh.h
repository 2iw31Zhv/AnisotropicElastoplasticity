#pragma once
#include <Eigen/Core>
#include <vector>

namespace igl
{
	namespace viewer
	{
		class Viewer;
	}
}

class LagrangianMesh
{
private:
	igl::viewer::Viewer * viewer_;
public:
	Eigen::MatrixX3d vertices;
	Eigen::MatrixX3d elements;

	Eigen::MatrixX3i faces;

	std::vector<Eigen::Matrix3d> plasticDeformationGradient;
	std::vector<Eigen::Matrix3d> elasticDeformationGradient;

	LagrangianMesh(const Eigen::MatrixX3d& vertices,
		const Eigen::MatrixX3i& faces,
		const std::vector<Eigen::Matrix3d>& elasticDeformationGradient,
		const std::vector<Eigen::Matrix3d>& plasticDeformationGradient);

	static LagrangianMesh ObjMesh(const std::string& filename);

	void bindViewer(igl::viewer::Viewer * viewer) { viewer_ = viewer; }
	void updateViewer();
};