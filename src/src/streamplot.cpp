/* Sourcefile for custom streamplotting  */

#include <cmath>
#include <vector>
#include <forward_list>
#include <list>
#include <array>
#include <exception>

#include <streamplot.hpp>

namespace acplotoo {

/* Vectorfield: */


/* Bilinear interpolation of the velocity grid at a point inside
 * the grid: */
std::array<double,2> interpolate(double x, double y, const Grid& velocity_grid,
                                 const Vectorfield& velocity)
{
	if (!velocity_grid.contains(x,y))
		throw std::runtime_error("Velocity grid does not contain point!");

	/* We know that velocity grid contains (x,y).
	 * Determine its four neighbours: */
	auto neighbours = velocity_grid.cell(x,y);

	/* Do bilinear velocity interpolation. */
	Point p0 = velocity_grid[neighbours[BOT_LEFT]];
	Point p1 = velocity_grid[neighbours[BOT_RIGHT]];
	Point p2 = velocity_grid[neighbours[TOP_LEFT]];
	Point p3 = velocity_grid[neighbours[TOP_RIGHT]];
	auto v0 = velocity[neighbours[BOT_LEFT]];
	auto v1 = velocity[neighbours[BOT_RIGHT]];
	auto v2 = velocity[neighbours[TOP_LEFT]];
	auto v3 = velocity[neighbours[TOP_RIGHT]];

	double cx = (x-p0.x()) / (p1.x() - p0.x());
	double cy = (y-p0.y()) / (p2.y() - p0.y());

	double vx1 = cx * v0.first + (1.0-cx) * v1.first;
	double vy1 = cx * v0.second + (1.0-cx) * v1.second;
	double vx2 = cx * v2.first + (1.0-cx) * v3.first;
	double vy2 = cx * v2.second + (1.0-cx) * v3.second;

	double vx = cy*vx1 + (1.0-cy)*vx2;
	double vy = cy*vy1 + (1.0-cy)*vy2;

	return {vx,vy};
}



/* Runge-Kutta integrator for generation of paths: */
std::vector<std::pair<double,double>>
_integrator_rk12(double x0, double y0, const Vectorfield& velocity,
                 const Grid& velocity_grid, Mask& mask,
                 const Grid& mask_grid, size_t id,
                 bool forward, bool backward, double tol,
                 double max_len, unsigned int max_steps=200)
{
	/* Initialize: */
	const std::array<bool,2> flag = {forward,backward};
	std::array<size_t,2> len = {0,0};
	std::array<std::forward_list<Point>,2> subtrajectories;
	double trajectory_length2 = 0;
	double max_len2 = max_len * max_len;

	for (int i : {0,1}){
		Grid::index_t current_tile = mask_grid.closest(x0, y0);
		auto last_tile = current_tile;
		if (!flag[i])
			continue;
		/* Do forward integration: */
		double x=x0;
		double y=y0;
		/* TODO Make this adjustable later on. */
		constexpr double dt0 = 1e-4;
		constexpr double dt_min = 1e-12;
		constexpr double dt_max = 1e-2;
		double dt = dt0;
		unsigned int step=0;
		while (velocity_grid.contains(x,y)
		       /* Also make sure that the current tile is not blocked by
		        * another trajectory (or that we have self-intersection): */
		       && (mask[current_tile] == 0
		           || (mask[current_tile] == id && current_tile == last_tile))
		       && step++ < max_steps
		       && trajectory_length2 < max_len2)
		{
			/* Set id: */
			mask[current_tile] = id;

			/******
			 * RKF45
			 * k1 = dt*f(t,x)
			 * k2 = dt*f(t+h/4, x+k/4)
			 * k3 = dt*f(t+3/8*h, x + 3/32*k1 + 9/32*k2)
			 * k4 = dt*f(t+12/13*h, x + 1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3)
			 * k5 = dt*f(t+h, x + 439/216*k1 - 8*k2 + 3680/513*k3 - 845/4104*k4)
			 * k6 = dt*f(t+h/2, x - 8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4
			 *           - 11/40*k5)
			 * Now we have a static vector field and can drop the time component.
			 * Follow:
			 * ece.uwaterloo.ca/~dwharder/NumericalAnalysis/14IVPs/rkf45/theory.html
			 */

			/* Obtain vector grid at current position: */
			auto v = interpolate(x, y, velocity_grid, velocity);
			double k1_x = dt * v[0];
			double k1_y = dt * v[1];

			/* k2: */
			v = interpolate(x + 0.5*k1_x, y + 0.5*k1_y, velocity_grid,
			                velocity);
			double k2_x = dt * v[0];
			double k2_y = dt * v[1];

			/* k3: */
			v = interpolate(x + 3./32.*k1_x + 9./32.*k2_x,
			                y + 3./32.*k1_y + 9./32.*k2_y,
			                velocity_grid, velocity);
			double k3_x = dt * v[0];
			double k3_y = dt * v[1];

			/* k4: */
			v = interpolate(x + 1932./2197.*k1_x - 7200./2197.*k2_x + 7296./2197.*k3_x,
			                y + 1932./2197.*k1_y - 7200./2197.*k2_y + 7296./2197.*k3_y,
			                velocity_grid, velocity);
			double k4_x = dt * v[0];
			double k4_y = dt * v[1];

			/* k5: */
			v = interpolate(x + 439./216.*k1_x - 8*k2_x + 3680./513.*k3_x
			                - 845./4104.*k4_x,
			                y + 439./216.*k1_y - 8*k2_y + 3680./513.*k3_y
			                - 845./4104.*k4_y,
			                velocity_grid, velocity);
			double k5_x = dt * v[0];
			double k5_y = dt * v[1];

			/* k6: */
			v = interpolate(x - 8./27.*k1_x + 2*k2_x - 3544./2565.*k3_x
			                + 1859./4104.*k4_x - 11./40.*k5_x,
			                y - 8./27.*k1_y + 2*k2_y - 3544./2565.*k3_y
			                + 1859./4104.*k4_y - 11./40.*k5_y,
			                velocity_grid, velocity);
			double k6_x = dt * v[0];
			double k6_y = dt * v[1];

			/* Final step for both RK4 and RK5: */
			double x_rk4 = x + 25./216.*k1_x + 1408./2565.*k3_x - 1./5.*k5_x;
			double y_rk4 = y + 25./216.*k1_y + 1408./2565.*k3_y - 1./5.*k5_y;

			double x_rk5 = x + 16./135.*k1_x + 6656./12825.*k3_x + 28561./56430.*k4_x
			               - 9./50.*k5_x + 2./55.*k6_x;
			double y_rk5 = y + 16./135.*k1_y + 6656./12825.*k3_y + 28561./56430.*k4_y
			               - 9./50.*k5_y + 2./55.*k6_y;

			/* Compare both results and calculate new step size: */
			double diff = x_rk4-x_rk5;
			double s = std::pow(0.5*tol*dt/std::abs(diff),0.25);
			if (s > 1.0){
				/* Accepted, use RK4 and adjust step size: */
				dt = std::min(s*dt,dt_max);

				/* Calculate new trajectory length: */
				trajectory_length2 += (x_rk4-x)*(x_rk4-x) + (y_rk4-y)*(y_rk4-y);

				/* Update coordinate: */
				x = x_rk4;
				y = y_rk4;
				last_tile = current_tile;
				current_tile = mask_grid.closest(x,y);

				/* Add to trajectory: */
				subtrajectories[i].emplace_front(x,y);
				++len[i];
			} else {
				/* Integration failed, reduce stepsize: */
				if (dt <= dt_min)
					/* If we're already at minimum step size, raise error: */
					throw std::runtime_error("Minimum step size did not converge!");
				dt = std::max(s*dt,dt_min);
			}
		}
	}

	/* Now we have the set of forward and backward trajectories. Save to vectors. */
	std::vector<std::pair<double,double>> trajectory(len[0]+len[1]+1);

	/* Between the subtrajectories, we save the starting point: */
	trajectory[len[0]] = {x0,y0};

	/* In reverse list order, add backward trajectory to front: */
	for (size_t i=0; i<len[0]; ++i){
		trajectory[len[0]-i-1] = subtrajectories[0].front();
		subtrajectories[0].pop_front();
	}
	/* Forward trajcectory to back: */
	for (size_t i=0; i<len[1]; ++i){
		trajectory[len[1]-i-1] = subtrajectories[1].front();
		subtrajectories[1].pop_front();
	}

	return trajectory;
}

std::list<std::vector<std::pair<double,double>>>
    streamplot_polygons(const std::vector<std::pair<double,double>>& start,
                        const Vectorfield& velocity, const Grid& velocity_grid,
                        const std::pair<size_t,size_t>& mask_size,
                        double min_len, double max_len, size_t max_steps,
                        bool forward, bool backward, double tol)
{
	/* Preparation: */
	const double min_len2 = min_len * min_len;
	/* Create the list of trajectories: */
	std::list<std::vector<std::pair<double,double>>> trajectories;

	/* Initialize mask! */
	Mask mask(mask_size.first, mask_size.second, 0);
	Grid mask_grid = velocity_grid.resample(mask_size.first, mask_size.second);

	for (size_t i=0; i<start.size(); ++i)
	{
		/* Integrate the trajectory: */
		std::vector<std::pair<double,double>> trajectory
		    = _integrator_rk12(start[i].first, start[i].second, velocity,
		                       velocity_grid, mask, mask_grid, i,
		                       forward, backward, tol, max_steps, max_len);

		/* If trajectory empty or too short, return: */
		if (trajectory.size() <= 2)
			continue;

		/* Calculate the trajectory distance and make sure it is long
		 * enough, otherwise continue with next trajectory: */
		double dist2 = 0.0;
		double x = trajectory[0].first;
		double y = trajectory[0].second;
		for (auto p : trajectory){
			dist2 += (x-p.first)*(x-p.first) + (y-p.second)*(y-p.second);
			x = p.first;
			y = p.second;
		}
		if (dist2 < min_len2)
			continue;

		/* If the trajectory is long enough, calculate the envelope: */
		/* TODO! */
		
		/* Now add the envelope to list: */
		trajectories.push_back(trajectory);
	}

	return trajectories;
}

} // End namespace

/* TODO remove! */
int main(){
	return 0;
}
