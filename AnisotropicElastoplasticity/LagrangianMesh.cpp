#include "LagrangianMesh.h"

#include <iostream>
#include <fstream>

#include <igl/viewer/Viewer.h>

#include "geometry.h"

using namespace std;
using namespace Eigen;

void LagrangianMesh::buildFaceWings_()
{
	int nfaces = faces.rows();
	faceWings.resize(nfaces, 3);
	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int result = -1;
			int p1 = faces(i, (j + 1) % 3);
			int p2 = faces(i, (j + 2) % 3);
			for (int k = 0; k < nfaces; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					if (faces(k, (l + 1) % 3) == p2 && faces(k, (l + 2) % 3) == p1)
					{
						result = faces(k, l);
					}
				}
			}
			faceWings(i, j) = result;
		}
	}
}

void LagrangianMesh::computeRestMetrics_()
{
	inverseMetrics.resize(faces.rows());

	for (int i = 0; i < faces.rows(); ++i)
	{
		Vector3d v1 = vertexPositions.row(faces(i, 0));
		Vector3d v2 = vertexPositions.row(faces(i, 1));
		Vector3d v3 = vertexPositions.row(faces(i, 2));

		double du1 = (v2 - v1).x();
		double du2 = (v3 - v1).x();

		double dv1 = (v2 - v1).y();
		double dv2 = (v3 - v1).y();

		Matrix2d metric;
		metric(0, 0) = du1;
		metric(0, 1) = du2;
		metric(1, 0) = dv1;
		metric(1, 1) = dv2;

		inverseMetrics[i] = metric.inverse();
	}
}

void LagrangianMesh::computeAreas_()
{
	areas.resize(faces.rows());
	areas.setZero();

	for (int i = 0; i < faces.rows(); ++i)
	{
		Vector3d v1 = vertexPositions.row(faces(i, 0));
		Vector3d v2 = vertexPositions.row(faces(i, 1));
		Vector3d v3 = vertexPositions.row(faces(i, 2));

		double area = 0.5 * ((v2 - v1).cross(v3 - v1)).norm();
		areas[i] = area;
	}
}

void LagrangianMesh::forceRelaxer(
	Eigen::Matrix2d & rotation, 
	Eigen::Matrix2d & symmetry, 
	double & J,
	double tolerance)
{
	if (fabs(symmetry(0, 0)) < tolerance)
	{
		symmetry(0, 0) = 0.0;
	}
	if (fabs(symmetry(1, 1)) < tolerance)
	{
		symmetry(1, 1) = 0.0;
	}
	if (fabs(symmetry(0, 1)) < tolerance)
	{
		symmetry(0, 1) = 0.0;
	}
	if (fabs(symmetry(1, 0)) < tolerance)
	{
		symmetry(1, 0) = 0.0;
	}

	if (fabs(rotation(0, 1)) < tolerance)
	{
		rotation(0, 1) = 0.0;
	}
	if (fabs(rotation(1, 0)) < tolerance)
	{
		rotation(1, 0) = 0.0;
	}
	if (fabs(rotation(1, 1) - 1.0) < tolerance)
	{
		rotation(1, 1) = 1.0;
	}
	if (fabs(rotation(0, 0) - 1.0) < tolerance)
	{
		rotation(0, 0) = 1.0;
	}

	if (fabs(J) < tolerance)
	{
		J = 0.0;
	}
}

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
	const Eigen::MatrixX3d& elementRestDirections_3,
	double mu,
	double lambda,
	double shearStiffness,
	double stiffness,
	double frictionCoeff) :
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
	elementRestDirections_1_(elementRestDirections_1),
	elementRestDirections_2_(elementRestDirections_2),
	elementRestDirections_3_(elementRestDirections_3),
	mu(mu),
	lambda(lambda),
	shearStiffness(shearStiffness),
	stiffness(stiffness),
	frictionCoeff(frictionCoeff)
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
	updateElementPositions();

	colors_.resize(faces.rows(), 3);
	colors_.setOnes();

	buildFaceWings_();
	computeAreas_();
	computeRestMetrics_();
}

LagrangianMesh LagrangianMesh::ObjMesh(const std::string & filename,
	double density,
	double thickness,
	double youngsModulus,
	double poissonRatio,
	double shearStiffness,
	double stiffness,
	double frictionAngleInDegree)
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

		elementRestDirections_1.row(f) = v2 - v1;
		elementRestDirections_2.row(f) = v3 - v1;
		elementRestDirections_3.row(f) = normal;
	}

	VectorXd vertexMasses = density * vertexVolumes;
	VectorXd elementMasses = density * elementVolumes;

	double lambda = youngsModulus * poissonRatio / (1.0 + poissonRatio) / (1.0 - 2.0 * poissonRatio);
	double mu = youngsModulus / 2.0 / (1.0 + poissonRatio);

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
		elementRestDirections_3,
		mu,
		lambda,
		shearStiffness,
		stiffness,
		tan(frictionAngleInDegree * igl::PI / 180.0));
}

void LagrangianMesh::updateViewer()
{
	viewer_->data.set_mesh(vertexPositions, faces);
	viewer_->data.set_colors(colors_);
}

void LagrangianMesh::updateElementPositions()
{
	for (int f = 0; f < faces.rows(); ++f)
	{
		const Vector3d& v1 = vertexPositions.row(faces.row(f)[0]);
		const Vector3d& v2 = vertexPositions.row(faces.row(f)[1]);
		const Vector3d& v3 = vertexPositions.row(faces.row(f)[2]);
		elementPositions.row(f) = (v1 + v2 + v3) / 3.0;
	}
}

void LagrangianMesh::computeVertexInPlaneForces(Eigen::MatrixX3d & vertexForces,
	std::vector<Eigen::Matrix2d>& inPlanePiolaKirhoffStresses)
{
	int Nv = vertexPositions.rows();
	vertexForces.resize(Nv, 3);
	vertexForces.setZero();

	int Ne = elementPositions.rows();
	inPlanePiolaKirhoffStresses.resize(Ne);

	for (int f = 0; f < Ne; ++f)
	{
		Matrix3d restDirectionMatrix;
		Vector3d D1, D2, D3;

		D1 = elementRestDirections_1().row(f);
		D2 = elementRestDirections_2().row(f);
		D3 = elementRestDirections_3().row(f);

		restDirectionMatrix.col(0) = D1;
		restDirectionMatrix.col(1) = D2;
		restDirectionMatrix.col(2) = D3;

		Matrix3d directionMatrix = elasticDeformationGradient[f] * restDirectionMatrix;

		Matrix3d Q, R;
		Vector3d q1, q2, q3;

		GramSchmidtOrthonomalization(Q, R, directionMatrix);

		q1 = Q.col(0);
		q2 = Q.col(1);
		q3 = Q.col(2);

		Matrix2d inPlaneR = R.block<2, 2>(0, 0);
		
		Matrix3d restQ, restR;
		GramSchmidtOrthonomalization(restQ, restR, restDirectionMatrix);

		Matrix2d restInPlaneR = restR.block<2, 2>(0, 0);
		Matrix2d invRest = inverseR(restInPlaneR);

		double i11 = invRest(0, 0),
			i12 = invRest(0, 1),
			i22 = invRest(1, 1);

		Matrix2d refInPlaneR = inPlaneR * invRest;

		Matrix2d invRefMulDet;
		invRefMulDet << R(1, 1), -R(0, 1),
			0.0, R(0, 0);

		JacobiSVD<Matrix2d> svd(refInPlaneR, ComputeFullU | ComputeFullV);
		Matrix2d rotation = svd.matrixU() * svd.matrixV().transpose();
		Matrix2d symmetry = svd.matrixV() * svd.singularValues().asDiagonal()
			* svd.matrixV().transpose() - Matrix2d::Identity();
		double J = refInPlaneR.diagonal().prod() - 1.0;
		forceRelaxer(rotation, symmetry, J, 1e-15);

		Matrix2d piolaKirhoffStress = 2.0 * mu * rotation * symmetry
			+ lambda * J * invRefMulDet.transpose();

		inPlanePiolaKirhoffStresses[f] = piolaKirhoffStress;

		double dr11 = piolaKirhoffStress(0, 0),
			dr12 = piolaKirhoffStress(0, 1),
			dr22 = piolaKirhoffStress(1, 1);

		Vector3d f2 = -(dr11 * i11 + dr12 * i12) * q1;
		Vector3d f3 = -dr12 * i22 * q1 - dr22 * i22 * q2;
		Vector3d f1 = -(f2 + f3);

		vertexForces.row(faces.row(f)[0]) += f1;
		vertexForces.row(faces.row(f)[1]) += f2;
		vertexForces.row(faces.row(f)[2]) += f3;
	}
}
