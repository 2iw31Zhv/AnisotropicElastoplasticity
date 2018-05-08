#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/viewer/Viewer.h>
#include "ParticleSystem.h"
#include "RegularGrid.h"
#include "LagrangianMesh.h"
#include "HybridSolver.h"
#include "LevelSet.h"
#include <thread>


using namespace std;
using namespace Eigen;
using namespace igl;

HybridSolver solver;

bool pre_draw(viewer::Viewer &viewer)
{
	viewer.data.clear();
	solver.updateViewer();
	return false;
}
void simulate()
{
	solver.solve(0.3, 100.0, 0.95);
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
	if (key == 's' || key == 'S')
	{
		thread t(simulate);
		t.detach();
	}

	return true;
}

int main()
{
	viewer::Viewer viewer;

	ParticleSystem ps = ParticleSystem::SandCylinder(
	Vector3d(0.0, 0.0, 0.02),
		0.5,
		2.0,
		1e4);

	double meshRes = 1;
	double gridLen = 1.0 / meshRes;

	//RegularGrid rg(Vector3d(- 1.2 * totalLen, - 0.4 * totalLen, -0.9 * totalLen),
	//	Vector3d(0.4 * totalLen, 0.4 * totalLen, 0.1 * totalLen),
	//	Vector3i(80, 40, 50));
	RegularGrid rg(Vector3d(- 9 * gridLen, - 3 * gridLen, - 9 * gridLen),
		Vector3d(3 * gridLen, 3 * gridLen, 3 * gridLen),
		Vector3i(12, 6, 12));

	VectorXd constraints;
	constraints.resize(9);
	constraints.setZero();
	constraints[0] = 1.0;
	constraints[2] = 1.0;

	LagrangianMesh mesh = LagrangianMesh::ObjMesh(
		"square_double_1.obj", 2e3, 1e-3, 200, 0.3, 0.0, 4e4, 0.0);

	mesh.bindConstraints(&constraints);

	//solver.setParticleSystem(&ps);
	solver.setRegularGrid(&rg);
	solver.setLagrangianMesh(&mesh);

	using namespace std::placeholders;
	LevelSet gls_m05 = bind(groundLevelSet, _1, -3.0);
	DLevelSet dgls_m05 = bind(DgroundLevelSet, _1, -3.0);
	//LevelSet w2g_110 = bind(wall2groundLevelSet, _1, 1.0, 1.0, 0.0);
	//DLevelSet dw2g_110 = bind(Dwall2groundLevelSet, _1, 1.0, 1.0, 0.0);
	solver.setLevelSet(gls_m05, dgls_m05);

	solver.bindViewer(&viewer);

	
	viewer.core.point_size = 3.0;
	viewer.core.background_color = Vector4f(0.0, 0.0, 0.0, 0.0);
	viewer.core.is_animating= true;
	viewer.core.show_lines = true;
	viewer.callback_key_down = &key_down;
	viewer.callback_pre_draw = &pre_draw;
	
	viewer.launch();
	return 0;
}