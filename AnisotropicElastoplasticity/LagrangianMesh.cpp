#include "LagrangianMesh.h"

#include <iostream>
#include <fstream>

#include <igl/viewer/Viewer.h>

#include "geometry.h"

using namespace std;
using namespace Eigen;

LagrangianMesh::LagrangianMesh(
	const Eigen::MatrixX3d & vertexPositions,
	const Eigen::MatrixX3i & faces,
	const Eigen::MatrixX3d & vertexVelocities,
	const Eigen::MatrixX3d & elementVelocities,
	const std::vector<Eigen::Matrix3d>& elasticDeformationGradient,
	const std::vector<Eigen::Matrix3d>& plasticDeformationGradient,
	const Eigen::VectorXd& vertexMasses,
	const Eigen::VectorXd& vertexVolumes,
	const Eigen::VectorXd& elementMasses,
	const Eigen::VectorXd& elementVolumes,
	const Eigen::MatrixX3d& elementRestDirections_1,
	const Eigen::MatrixX3d& elementRestDirections_2,
	const Eigen::MatrixX3d& elementRestDirections_3) :
	vertexPositions(vertexPositions),
	faces(faces),
	vertexVelocities(vertexVelocities),
	elementVelocities(elementVelocities),
	elasticDeformationGradient(elasticDeformationGradient),
	plasticDeformationGradient(plasticDeformationGradient),
	vertexMasses(vertexMasses),
	vertexVolumes(vertexVolumes),
	elementMasses(elementMasses),
	elementVolumes(elementVolumes),
	elementRestDirections_1(elementRestDirections_1),
	elementRestDirections_2(elementRestDirections_2),
	elementRestDirections_3(elementRestDirections_3)
{
	int Nv = vertexPositions.rows();
	int Nf = faces.rows();

	elementPositions.resize(Nf, 3);
	elementPositions.setZero();

	vertexAffineMomenta_1.resize(Nv, 3);
	vertexAffineMomenta_2.resize(Nv, 3);
	vertexAffineMomenta_3.resize(Nv, 3);

	elementAffineMomenta_1.resize(Nf, 3);
	elementAffineMomenta_2.resize(Nf, 3);
	elementAffineMomenta_3.resize(Nf, 3);

	// set the element positions from the vertices
	for (int f = 0; f < Nf; ++f)
	{
		const Vector3d& v1 = vertexPositions.row(faces.row(f)[0]);
		const Vector3d& v2 = vertexPositions.row(faces.row(f)[1]);
		const Vector3d& v3 = vertexPositions.row(faces.row(f)[2]);
		elementPositions.row(f) = (v1 + v2 + v3) / 3.0;
	}

	colors_.resize(faces.rows(), 3);
	colors_.setOnes();


}

LagrangianMesh LagrangianMesh::ObjMesh(const std::string & filename,
	double density,
	double thickness)
{
	ifstream fin(filename);

	if (!fin.is_open())
	{
		cerr << "[ERROR]: Cannot open the obj file "
			<< filename << "! " << __FUNCTION__ << ", "
			<< __LINE__ << endl;
		exit(1);
	}

	std::string mark;
	float position[3];
	unsigned short indices[3];
	char buffer[1024];

	vector<double> posBuf;
	vector<int> faceBuf;
	// currently, we only support triangulated facets
	
	while (fin >> mark)
	{
		if (mark == "v")
		{
			fin >> position[0] >> position[1] >> position[2];
			posBuf.push_back(position[0]);
			posBuf.push_back(position[1]);
			posBuf.push_back(position[2]);

		}
		else if (mark == "vn")
		{
			fin.getline(buffer, 1024);
		}
		else if (mark == "vt")
		{
			fin.getline(buffer, 1024);
		}
		else if (mark == "f")
		{
			fin >> indices[0] >> indices[1] >> indices[2];
			faceBuf.push_back(indices[0]);
			faceBuf.push_back(indices[1]);
			faceBuf.push_back(indices[2]);

		}
		else if (mark == "#")
		{
			fin.getline(buffer, 1024);
		}
		else
		{
			fin.getline(buffer, 1024);
		}
	}
	
	MatrixX3i F;
	MatrixX3d V;
	F.resize(faceBuf.size() / 3, 3);
	V.resize(posBuf.size() / 3, 3);

	int fid = 0;
	int vid = 0;
	for (int f = 0; f < F.rows(); ++f)
	{
		F(f, 0) = faceBuf[fid++] - 1;
		F(f, 1) = faceBuf[fid++] - 1;
		F(f, 2) = faceBuf[fid++] - 1;
	}

	for (int i = 0; i < V.rows(); ++i)
	{
		V(i, 0) = posBuf[vid++];
		V(i, 1) = posBuf[vid++];
		V(i, 2) = posBuf[vid++];
	}

	MatrixX3d vertexVelocities;
	vertexVelocities.resize(V.rows(), 3);
	vertexVelocities.setZero();

	MatrixX3d elementVelocities;
	elementVelocities.resize(F.rows(), 3);
	elementVelocities.setZero();

	vector<Matrix3d> elasticDeformationGradient;
	vector<Matrix3d> plasticDeformationGradient;


	VectorXd vertexVolumes;
	vertexVolumes.resize(V.rows());
	vertexVolumes.setZero();

	VectorXd elementVolumes;
	elementVolumes.resize(F.rows());
	elementVolumes.setZero();

	MatrixX3d elementRestDirections_1;
	elementRestDirections_1.resize(F.rows(), 3);
	MatrixX3d elementRestDirections_2;
	elementRestDirections_2.resize(F.rows(), 3);
	MatrixX3d elementRestDirections_3;
	elementRestDirections_3.resize(F.rows(), 3);

	for (int f = 0; f < F.rows(); ++f)
	{
		elasticDeformationGradient.push_back(Matrix3d::Identity());
		plasticDeformationGradient.push_back(Matrix3d::Identity());

		const Vector3d& v1 = V.row(F.row(f)[0]);
		const Vector3d& v2 = V.row(F.row(f)[1]);
		const Vector3d& v3 = V.row(F.row(f)[2]);

		double area = triangleArea(v1, v2, v3);

		double volume = 0.25 * area * thickness;
		elementVolumes[f] = volume;

		vertexVolumes[F.row(f)[0]] += volume;
		vertexVolumes[F.row(f)[1]] += volume;
		vertexVolumes[F.row(f)[2]] += volume;

		Vector3d normal = triangleNormal(v1, v2, v3);
		Vector3d t1 = (v2 - v1).normalized();
		Vector3d t2 = normal.cross(t1);

		elementRestDirections_1.row(f) = t1;
		elementRestDirections_2.row(f) = t2;
		elementRestDirections_3.row(f) = normal;

	}

	VectorXd vertexMasses = density * vertexVolumes;
	VectorXd elementMasses = density * elementVolumes;

	return LagrangianMesh(V, F, 
		vertexVelocities,
		elementVelocities,
		elasticDeformationGradient, 
		plasticDeformationGradient,
		vertexMasses,
		vertexVolumes,
		elementMasses,
		elementVolumes,
		elementRestDirections_1,
		elementRestDirections_2,
		elementRestDirections_3);
}

void LagrangianMesh::updateViewer()
{
	viewer_->data.set_mesh(vertexPositions, faces);
	viewer_->data.set_colors(colors_);
}
