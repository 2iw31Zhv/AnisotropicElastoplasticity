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

	double meshRes = 28.4;
	double gridLen = 1.0 / meshRes;
	double totalLen = 50.0 * gridLen;

	RegularGrid rg(Vector3d(- totalLen, -totalLen, - totalLen),
		Vector3d(totalLen, totalLen, 3.0 * totalLen),
		Vector3i(100, 100, 200));
	
	LagrangianMesh mesh = LagrangianMesh::ObjMesh(
		"square_random.obj", 2e3, 1e-3, 200, 0.3, 0.0, 4e4, 0.0);

	solver.setParticleSystem(&ps);
	solver.setRegularGrid(&rg);
	solver.setLagrangianMesh(&mesh);

	using namespace std::placeholders;
	LevelSet gls_m05 = bind(groundLevelSet, _1, -0.5);
	DLevelSet dgls_m05 = bind(DgroundLevelSet, _1, -0.5);
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