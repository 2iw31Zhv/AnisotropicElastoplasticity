#include "LevelSet.h"


double groundLevelSet(const Eigen::Vector3d& x, double groundZ)
{
	return x[2] - groundZ;
}

Eigen::Vector3d DgroundLevelSet(const Eigen::Vector3d & x, double groundZ)
{
	return Eigen::Vector3d(0.0, 0.0, 1.0);
}
