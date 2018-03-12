#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/viewer/Viewer.h>
#include "ParticleSystem.h"
#include "RegularGrid.h"

using namespace std;
using namespace Eigen;
using namespace igl;

ParticleSystem ps = ParticleSystem::Ball(Vector3d(0.0, 0.0, 0.0), 1.0, 1000);
RegularGrid rg(Vector3d(-2.0, -2.0, -2.0),
	Vector3d(2.0, 2.0, 2.0),
	Vector3i(50, 50, 50));

bool pre_draw(viewer::Viewer &viewer)
{
	viewer.data.clear();
	ps.updateViewer();
	rg.updateViewer();
	return false;
}

int main()
{
	viewer::Viewer viewer;

	ps.bindViewer(&viewer);
	rg.bindViewer(&viewer);
	viewer.core.background_color = Vector4f(0.0, 0.0, 0.0, 0.0);
	viewer.core.is_animating= true;
	viewer.core.show_lines = true;
	//viewer.callback_key_down = &key_down;
	viewer.callback_pre_draw = &pre_draw;
	viewer.launch();
	return 0;
}