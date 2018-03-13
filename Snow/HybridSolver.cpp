#include "HybridSolver.h"
#include "ParticleSystem.h"
#include "RegularGrid.h"
#include <vector>
#include <igl/viewer/Viewer.h>
#include "interpolation.h"

using namespace std;
using namespace Eigen;
using namespace igl;

void HybridSolver::evaluateInterpolationWeights_()
{
	vector<Triplet<double>> omegaTriplets;
	for (int p = 0; p < ps_->positions.rows(); ++p)
	{
		const Vector3d& position = ps_->positions.row(p);
		const Vector3d& posRef = position - rg_->minBound();
		int i_fl = static_cast<int>(posRef[0] / rg_->h()[0]);
		int j_fl = static_cast<int>(posRef[1] / rg_->h()[1]);
		int k_fl = static_cast<int>(posRef[2] / rg_->h()[2]);

		for (int i = i_fl - 2; i <= i_fl + 2; ++i)
		{
			for (int j = j_fl - 2; j <= j_fl + 2; ++j)
			{
				for (int k = k_fl - 2; k <= k_fl + 2; ++k)
				{
					if (0 <= i && i < rg_->resolution()[0] &&
						0 <= j && j < rg_->resolution()[1] &&
						0 <= k && k < rg_->resolution()[2])
					{
						double dyadic_i = cubic_B_spline(posRef[0] / rg_->h()[0] - i);
						double dyadic_j = cubic_B_spline(posRef[1] / rg_->h()[1] - j);
						double dyadic_k = cubic_B_spline(posRef[2] / rg_->h()[2] - k);

						if (dyadic_i > 0 && dyadic_j > 0 && dyadic_k > 0)
						{
							omegaTriplets.push_back(Triplet<double>(p,
								rg_->toIndex(i, j, k),
								dyadic_i * dyadic_j * dyadic_k));
						}
					}
				}
			}
		}
	}

	omegas_.setZero();
	omegas_.setFromTriplets(omegaTriplets.begin(), omegaTriplets.end());
}

void HybridSolver::solve(double Dt)
{
	// rasterize particle data to grid: m, v

	// compute particle densities and volumes

	// while true
	//	compute grid force

	//	update velocities

	//	collision projection

	//	update deformation gradient
	
	//	update particle velocities

	//	update particle positions

	//	rasterize particle data to grid: m, v
}

void HybridSolver::bindViewer(igl::viewer::Viewer * viewer)
{
	viewer_ = viewer;
	if (ps_ != nullptr)
	{
		ps_->bindViewer(viewer);
	}

	if (rg_ != nullptr)
	{
		rg_->bindViewer(viewer);
	}
}

void HybridSolver::updateViewer()
{
	using namespace viewer;
	ps_->updateViewer();
	rg_->updateViewer();
}
