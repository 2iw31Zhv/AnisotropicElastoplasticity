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
	Vector3d(-0.5, 0.0, 0.1),
		0.25,
		0.6,
		1e4);

	double meshRes = 24.0;
	double gridLen = 1.0 / meshRes;
	double totalLen = 50.0 * gridLen;

	RegularGrid rg(Vector3d(-1.2 * totalLen, - 0.6 * totalLen, -0.8 * totalLen),
		Vector3d(0.6 * totalLen, 0.6 * totalLen, 0.8 * totalLen),
		Vector3i(90, 60, 80));

	VectorXd constraints;
	constraints.resize(425);
	constraints.setZero();
	constraints[0] = 1.0;
	constraints[16] = 1.0;

	LagrangianMesh mesh = LagrangianMesh::ObjMesh(
		"square_double_12x8.obj", 2e3, 0.04, 200, 0.3, 0.0, 4e4, 0.0);

	mesh.bindConstraints(&constraints);

	solver.setParticleSystem(&ps);
	solver.setRegularGrid(&rg);
	solver.setLagrangianMesh(&mesh);

	using namespace std::placeholders;
	LevelSet gls_m14 = bind(groundLevelSet, _1, -1.4);
	DLevelSet dgls_m14 = bind(DgroundLevelSet, _1, -1.4);
	//LevelSet w2g_110 = bind(wall2groundLevelSet, _1, 1.0, 1.0, 0.0);
	//DLevelSet dw2g_110 = bind(Dwall2groundLevelSet, _1, 1.0, 1.0, 0.0);
	solver.setLevelSet(gls_m14, dgls_m14);

	solver.bindViewer(&viewer);

	
	viewer.core.point_size = 3.0;
	viewer.core.background_color = Vector4f(0.0, 0.0, 0.0, 0.0);
	viewer.core.is_animating= true;
	viewer.core.show_lines = false;
	viewer.callback_key_down = &key_down;
	viewer.callback_pre_draw = &pre_draw;
	
	viewer.launch();
	return 0;
}