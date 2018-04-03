#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/viewer/Viewer.h>
#include "ParticleSystem.h"
#include "RegularGrid.h"
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
	solver.solve(1e-3, 100.0, 0.95);
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

	ParticleSystem ps = ParticleSystem::SandBlock(
	Vector3d(0.0, 0.0, 0.0),
		Vector3d(1.0, 1.0, 2.0),
		0.5,
		1e4);

	RegularGrid rg(Vector3d(-7.5, -7.5, -0.5),
		Vector3d(2.5, 2.5, 7.5),
		Vector3i(100, 100, 80));

	solver.setParticleSystem(&ps);
	solver.setRegularGrid(&rg);
	

	using namespace std::placeholders;
	//LevelSet gls_m00 = bind(groundLevelSet, _1, 0.0);
	//DLevelSet dgls_m00 = bind(DgroundLevelSet, _1, 0.0);
	LevelSet w2g_110 = bind(wall2groundLevelSet, _1, 1.0, 1.0, 0.0);
	DLevelSet dw2g_110 = bind(Dwall2groundLevelSet, _1, 1.0, 1.0, 0.0);
	solver.setLevelSet(w2g_110, dw2g_110);

	solver.bindViewer(&viewer);


	viewer.core.point_size = 1.0;
	viewer.core.background_color = Vector4f(0.0, 0.0, 0.0, 0.0);
	viewer.core.is_animating= true;
	viewer.core.show_lines = true;
	viewer.core.light_position = Vector3f(-0.3, -0.3, 3.0);
	viewer.callback_key_down = &key_down;
	viewer.callback_pre_draw = &pre_draw;
	viewer.launch();
	return 0;
}