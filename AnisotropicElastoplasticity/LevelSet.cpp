#include "LevelSet.h"

#include <algorithm>

using namespace std;
using namespace Eigen;

double groundLevelSet(const Eigen::Vector3d& x, double groundZ)
{
	return x[2] - groundZ;
}

Eigen::Vector3d DgroundLevelSet(const Eigen::Vector3d & x, double groundZ)
{
	return Eigen::Vector3d(0.0, 0.0, 1.0);
}

double wall2groundLevelSet(const Eigen::Vector3d & x, double wallX, double wallY, double groundZ)
{
	return min(min(x[2] - groundZ, wallX - x[0]), wallY - x[1]);
}

Eigen::Vector3d Dwall2groundLevelSet(const Eigen::Vector3d & x, double wallX, double wallY, double groundZ)
{
	// check the closest wall
	double dz = abs(x[2] - groundZ);
	double dx = abs(wallX - x[0]);
	double dy = abs(wallY - x[1]);

	if (dz <= dx && dz <= dy)
	{
		return Vector3d(0.0, 0.0, 1.0);
	}
	else if (dy <= dx)
	{
		return Vector3d(0.0, -1.0, 0.0);
	}
	else
	{
		return Vector3d(-1.0, 0.0, 0.0);
	}
}
