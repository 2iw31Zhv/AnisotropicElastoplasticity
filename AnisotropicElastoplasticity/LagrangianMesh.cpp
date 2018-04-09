#include "LagrangianMesh.h"

#include <iostream>
#include <fstream>

#include <igl/viewer/Viewer.h>

using namespace std;
using namespace Eigen;

LagrangianMesh::LagrangianMesh(
	const Eigen::MatrixX3d & vertices, 
	const Eigen::MatrixX3i & faces, 
	const std::vector<Eigen::Matrix3d>& elasticDeformationGradient, 
	const std::vector<Eigen::Matrix3d>& plasticDeformationGradient):
	vertices(vertices),
	faces(faces),
	elasticDeformationGradient(elasticDeformationGradient),
	plasticDeformationGradient(plasticDeformationGradient)
{
	// set the element positions from the vertices
	elements.resize(faces.rows(), 3);
	elements.setZero();
	for (int f = 0; f < faces.rows(); ++f)
	{
		const Vector3d& v1 = vertices.row(faces.row(f)[0]);
		const Vector3d& v2 = vertices.row(faces.row(f)[1]);
		const Vector3d& v3 = vertices.row(faces.row(f)[2]);
		elements.row(f) = (v1 + v2 + v3) / 3.0;
	}
}

LagrangianMesh LagrangianMesh::ObjMesh(const std::string & filename)
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

	vector<Matrix3d> elasticDeformationGradient;
	vector<Matrix3d> plasticDeformationGradient;
	for (int f = 0; f < F.rows(); ++f)
	{
		elasticDeformationGradient.push_back(Matrix3d::Identity());
		plasticDeformationGradient.push_back(Matrix3d::Identity());
	}
	return LagrangianMesh(V, F, elasticDeformationGradient, plasticDeformationGradient);
}

void LagrangianMesh::updateViewer()
{
	viewer_->data.set_mesh(vertices, faces);
}
