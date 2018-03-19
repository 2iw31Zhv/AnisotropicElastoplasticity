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
	solver.solve(0.001, 100.0, 0.95);
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

	ParticleSystem ps = ParticleSystem::SnowBall(Vector3d(0.0, 0.0, 0.0), 1.0, 1000);
	RegularGrid rg(Vector3d(-2.0, -2.0, -2.0),
		Vector3d(2.0, 2.0, 2.0),
		Vector3i(50, 50, 50));

	solver.setParticleSystem(&ps);
	solver.setRegularGrid(&rg);
	

	using namespace std::placeholders;
	LevelSet gls_m18 = bind(groundLevelSet, _1, -1.8);
	DLevelSet dgls_m18 = bind(DgroundLevelSet, _1, -1.8);
	solver.setLevelSet(gls_m18, dgls_m18);

	solver.bindViewer(&viewer);
	
	viewer.core.background_color = Vector4f(0.0, 0.0, 0.0, 0.0);
	viewer.core.is_animating= true;
	viewer.core.show_lines = true;
	
	viewer.callback_key_down = &key_down;
	viewer.callback_pre_draw = &pre_draw;
	viewer.launch();
	return 0;
}