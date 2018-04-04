#pragma once
#include <Eigen/Core>
#include <functional>

using LevelSet = std::function<double(const Eigen::Vector3d&)>;
using DLevelSet = std::function<Eigen::Vector3d(const Eigen::Vector3d&)>;

double groundLevelSet(const Eigen::Vector3d& x, double groundZ);
Eigen::Vector3d DgroundLevelSet(const Eigen::Vector3d& x, double groundZ);

double wall2groundLevelSet(const Eigen::Vector3d& x, double wallX, double wallY, double groundZ);
Eigen::Vector3d Dwall2groundLevelSet(const Eigen::Vector3d& x, double wallX, double wallY, double groundZ);