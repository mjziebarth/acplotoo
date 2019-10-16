/* Sourcefile for custom streamplotting  */

#include <cmath>
#include <vector>
#include <forward_list>
#include <list>
#include <array>
#include <exception>
#include <deque>

#include <streamplot.hpp>
#include <pointset.hpp>

namespace acplotoo {

/* Bilinear interpolation of the velocity grid at a point inside
 * the grid: */
Rect<Grid::index_t> interpolation_coefficients(double x, double y, const Grid& grid, 
                                               double& cx, double& cy)
{
	if (!grid.contains(x,y))
		throw std::runtime_error("Velocity grid does not contain point!");

	/* We know that velocity grid contains (x,y).
	 * Determine its four neighbours: */
	auto neighbours = grid.cell(x,y);

	/* Do bilinear velocity interpolation. */
	Point p0 = grid[neighbours[BOT_LEFT]];
	Point p1 = grid[neighbours[BOT_RIGHT]];
	Point p2 = grid[neighbours[TOP_LEFT]];

	if (p1.x() == p0.x()){
		cx = 0.5;
	} else {
		cx = (x-p0.x()) / (p1.x() - p0.x());
	}
	if (p1.y() == p0.y()){
		cy = 0.5;
	} else {
		cy = (y-p0.y()) / (p2.y() - p0.y());
	}

	return neighbours;
}


std::array<double,2> interpolate(double x, double y, const Grid& velocity_grid,
                                 const Vectorfield& velocity)
{
	double cx, cy;
	Rect<Grid::index_t> neighbours = interpolation_coefficients(x, y, velocity_grid,
	                                                            cx, cy);

	auto v0 = velocity[neighbours[BOT_LEFT]];
	auto v1 = velocity[neighbours[BOT_RIGHT]];
	auto v2 = velocity[neighbours[TOP_LEFT]];
	auto v3 = velocity[neighbours[TOP_RIGHT]];

	double vx1 = cx * v0.first + (1.0-cx) * v1.first;
	double vy1 = cx * v0.second + (1.0-cx) * v1.second;
	double vx2 = cx * v2.first + (1.0-cx) * v3.first;
	double vy2 = cx * v2.second + (1.0-cx) * v3.second;

	double vx = cy*vx1 + (1.0-cy)*vx2;
	double vy = cy*vy1 + (1.0-cy)*vy2;

	return {vx,vy};
}

double interpolate(double x, double y, const Grid& grid, const Scalarfield& scalar)
{
	double cx, cy;
	Rect<Grid::index_t> neighbours = interpolation_coefficients(x, y, grid, cx, cy);

	double s0 = scalar[neighbours[BOT_LEFT]];
	double s1 = scalar[neighbours[BOT_RIGHT]];
	double s2 = scalar[neighbours[TOP_LEFT]];
	double s3 = scalar[neighbours[TOP_RIGHT]];

	double s01 = cx * s0 + (1.0-cx) * s1;
	double s23 = cx * s2 + (1.0-cx) * s3;

	return cy*s01 + (1.0-cy)*s23;
}


bool check_collision(double x, double y, double r, size_t id,
                     const PointSet& points)
{
	/* Obtain all within range of (x,y) and see whether there are any
	 * of other ID: */
	auto in_range = points.query_in_range(x, y, r);
	for (auto& pt : in_range){
		/* Obtain id: */
		unsigned short id_i = std::get<2>(pt);

		/* See if the element id indicates a collision: */
		if (id_i && id_i != id){
			return false;
		}
	}

	return true;
}



double vector_length(const std::pair<double,double>& p)
{
	return std::sqrt(p.first*p.first + p.second*p.second);
}

/* Runge-Kutta integrator for generation of paths: */
std::pair<std::vector<std::pair<double,double>>,double>
_integrator_rk45(double x0, double y0, const Vectorfield& velocity,
                 const Grid& velocity_grid,
                 const PointSet& points, const tid_t id,
                 bool forward, bool backward, double tol,
                 double min_len, double max_len, double dt_min, double dt_max,
                 double step_len_min, double collision_radius,
                 unsigned int max_steps, unsigned int tile_history_size)
{
	/* Initialize: */
	const std::array<bool,2> flag = {forward,backward};
	std::array<size_t,2> len = {0,0};
	std::array<std::forward_list<Point>,2> subtrajectories;
	double trajectory_length = 0.0;
	double traj_len_at_last_tick = 0.0;
	bool leaves_rect = false;

	for (int i : {0,1}){
		/* Decide whether we should integrate in the current direction: */
		if (!flag[i])
			continue;
		/* Do forward integration: */
		double x=x0;
		double y=y0;
		/* TODO Make this adjustable later on. */
		const double dt0 = 0.1 * tol / vector_length(velocity[
		                                     velocity_grid.closest(x0,y0)]);
		double dt = dt0;
		double dt_sign = (i == 0) ? 1 : -1;
		unsigned int step=0;
		bool recalibrate = false;
		while (velocity_grid.contains(x,y)
		       /* Also make sure that the current tile is not blocked by
		        * another trajectory (or that we have self-intersection): */
		       && (recalibrate ||
		           check_collision(x, y, collision_radius, id, points))
		       && step++ < max_steps
		       && trajectory_length < max_len)
		{

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
			double k1_x = dt * v[0] * dt_sign;
			double k1_y = dt * v[1] * dt_sign;

			/* k2: */
			double x_test = x + 0.25*k1_x;
			double y_test = y + 0.25*k1_y;
			if (!velocity_grid.contains(x_test,y_test))
				continue;
			v = interpolate(x_test, y_test, velocity_grid,
			                velocity);
			double k2_x = dt * v[0] * dt_sign;
			double k2_y = dt * v[1] * dt_sign;

			/* k3: */
			x_test = x + 3./32.*k1_x + 9./32.*k2_x;
			y_test = y + 3./32.*k1_y + 9./32.*k2_y;
			if (!velocity_grid.contains(x_test,y_test))
				continue;
			v = interpolate(x_test, y_test, velocity_grid, velocity);
			double k3_x = dt * v[0] * dt_sign;
			double k3_y = dt * v[1] * dt_sign;

			/* k4: */
			x_test = x + 1932./2197.*k1_x - 7200./2197.*k2_x + 7296./2197.*k3_x;
			y_test = y + 1932./2197.*k1_y - 7200./2197.*k2_y + 7296./2197.*k3_y;
			if (!velocity_grid.contains(x_test,y_test))
				continue;
			v = interpolate(x_test, y_test, velocity_grid, velocity);
			double k4_x = dt * v[0] * dt_sign;
			double k4_y = dt * v[1] * dt_sign;

			/* k5: */
			x_test = x + 439./216.*k1_x - 8*k2_x + 3680./513.*k3_x
			         - 845./4104.*k4_x;
			y_test = y + 439./216.*k1_y - 8*k2_y + 3680./513.*k3_y
			         - 845./4104.*k4_y;
			if (!velocity_grid.contains(x_test,y_test))
				continue;
			v = interpolate(x_test, y_test, velocity_grid, velocity);
			double k5_x = dt * v[0] * dt_sign;
			double k5_y = dt * v[1] * dt_sign;

			/* k6: */
			x_test = x - 8./27.*k1_x + 2*k2_x - 3544./2565.*k3_x
			                + 1859./4104.*k4_x - 11./40.*k5_x;
			y_test = y - 8./27.*k1_y + 2*k2_y - 3544./2565.*k3_y
			                + 1859./4104.*k4_y - 11./40.*k5_y;
			if (!velocity_grid.contains(x_test,y_test))
				continue;
			v = interpolate(x_test, y_test, velocity_grid, velocity);
			double k6_x = dt * v[0] * dt_sign;
			double k6_y = dt * v[1] * dt_sign;

			/* Final step for both RK4 and RK5: */
			double x_rk4 = x + 25./216.*k1_x + 1408./2565.*k3_x + 2197./4101.*k4_x
			                 - 1./5.*k5_x;
			double y_rk4 = y + 25./216.*k1_y + 1408./2565.*k3_y + 2197./4101.*k4_y
			                 - 1./5.*k5_y;

			double x_rk5 = x + 16./135.*k1_x + 6656./12825.*k3_x + 28561./56430.*k4_x
			               - 9./50.*k5_x + 2./55.*k6_x;
			double y_rk5 = y + 16./135.*k1_y + 6656./12825.*k3_y + 28561./56430.*k4_y
			               - 9./50.*k5_y + 2./55.*k6_y;

			/* Compare both results and calculate new step size: */
			double diff = std::sqrt((x_rk4-x_rk5)*(x_rk4-x_rk5) + 
			                        (y_rk4-y_rk5)*(y_rk4-y_rk5));
			double s = std::pow(0.5*tol*dt/diff,0.25);
			if (s > 1.0){
				/* Accepted, use RK4 and adjust step size: */
				dt = std::min(s*dt,dt_max);
				recalibrate = false; // Actually only measures whether non-accepted. 

				/* Calculate new trajectory length: */
				trajectory_length += std::sqrt((x_rk4-x)*(x_rk4-x)
				                               + (y_rk4-y)*(y_rk4-y));

				/* Update coordinate: */
				x = x_rk4;
				y = y_rk4;

				/* Add to trajectory, if far enough away:: */
				if (trajectory_length - traj_len_at_last_tick >= step_len_min){
					traj_len_at_last_tick = trajectory_length;
					subtrajectories[i].emplace_front(x,y);
					++len[i];
				}
			} else {
				/* Integration failed, reduce stepsize: */
				if (dt <= dt_min)
					/* If we're already at minimum step size, raise error: */
					throw std::runtime_error("Minimum step size did not converge!");
				dt = std::max(s*dt,dt_min);
				recalibrate = true;
			}
		}
		leaves_rect = leaves_rect || !velocity_grid.contains(x,y);
	}

	/* If trajectory length too small, return empty trajectory: */
	if (!leaves_rect && trajectory_length < min_len){
		std::pair<std::vector<std::pair<double,double>>,double> retval;
		retval.second = 0;
		return retval;
	}

	/* Now we have the set of forward and backward trajectories. Save to vectors. */
	std::pair<std::vector<std::pair<double,double>>,double> retval;
	std::vector<std::pair<double,double>>& trajectory(retval.first);
	trajectory.resize(len[0]+len[1]+1);
	
	/* Save length: */
	retval.second = trajectory_length;

	/* Between the subtrajectories, we save the starting point: */
	trajectory[len[1]] = {x0,y0};

	/* In reverse list order, add backward trajectory to front.
	 * At the same time, add to points. */
	for (size_t i=0; i<len[1]; ++i){
		trajectory[i] = subtrajectories[1].front();
		subtrajectories[1].pop_front();
	}

	/* Forward trajectory to back: */
	for (size_t i=0; i<len[0]; ++i){
		trajectory[len[0]+len[1]-i] = subtrajectories[0].front();
		subtrajectories[0].pop_front();
	}

	return retval;
}

std::pair<std::list<std::vector<std::pair<double,double>>>,
          std::vector<std::array<double,4>>>
    streamplot_polygons(const std::vector<std::pair<double,double>>& start,
                        const Vectorfield& velocity, const Grid& velocity_grid,
                        const Scalarfield& width_field,
                        double min_len, double max_len, double step_len_min,
                        double arrow_head_step, double collision_radius,
                        size_t max_steps, bool forward, bool backward, double tol,
                        size_t tile_history_size)
{
	/* Make sure that number of starting points does not exceed tid_t range: */
	if (start.size()+1 >= static_cast<tid_t>(-1))
	{
		throw std::runtime_error("Too many starting points for chosen "
		                         "trajectory id data type.");
	}

	/* Create the list of trajectories: */
	std::list<std::vector<std::pair<double,double>>> trajectories;
	std::vector<std::array<double,4>> all_heads;

	/* Determine optimal step sizes: */
	double max_v = 0.0;
	for (size_t i=0; i<velocity.nx(); ++i){
		for (size_t j=0; j<velocity.ny(); ++j){
			auto p = velocity[{i,j}];
			double v = p.first * p.first + p.second * p.second;
			if (v > max_v)
				max_v = v;
		}
	}
	max_v = std::sqrt(max_v);
	double dt_max = 0.01 / max_v;
	double dt_min = 1e-40 / max_v;


	/* Initialize the set of trajectory points: */
	PointSet points;

	for (size_t i=0; i<start.size(); ++i)
	{
		/* Integrate the trajectory: */
		auto retval
		    = _integrator_rk45(start[i].first, start[i].second, velocity,
		                       velocity_grid, points, i+1,
		                       forward, backward, tol, min_len, max_len,
		                       dt_min, dt_max, step_len_min, collision_radius,
		                       max_steps, tile_history_size);

		std::vector<std::pair<double,double>>& trajectory = retval.first;
		const double whole_trajectory_len = retval.second;

		/* If trajectory empty or too short, return: */
		if (trajectory.size() <= 2)
			continue;

		/* Determine the widths along the trajectory. */
		std::vector<double> widths(trajectory.size());
		for (size_t j=0; j<trajectory.size(); ++j){
			if (!velocity_grid.contains(trajectory[j].first,
			                            trajectory[j].second))
			{
				widths[j] = std::nan("");
			} else {
				widths[j] = interpolate(trajectory[j].first, trajectory[j].second,
				                        velocity_grid, width_field);
			}
		}

		/* If widths are NaN, we split the trajectory. */
		struct trajectory_t {
			size_t j0;
			size_t j1;
			trajectory_t(size_t j0, size_t j1) : j0(j0),j1(j1) {};
		};
		std::vector<trajectory_t> subtrajectories;
		size_t j0 = 0;
		for (size_t j=0; j<trajectory.size(); ++j){
			if (std::isnan(widths[j])){
				/* Split the trajectory and add subtrajectory
				 * if long enough: */
				if (j - j0 > 2){
					subtrajectories.emplace_back(j0, j);
				}
				j0 = j;
			}
		}
		if (trajectory.size() - j0 > 2)
			subtrajectories.emplace_back(j0,trajectory.size());

		if (subtrajectories.empty())
			continue;

		for (auto t : subtrajectories){
			size_t len = t.j1 - t.j0;
			/* Determine the arrow head positions: */
			double trajectory_len = (whole_trajectory_len * trajectory.size())
			                        / len;
			size_t n_heads = std::max<size_t>(std::floor(trajectory_len /
			                                             arrow_head_step)+1,
			                                  2) - 1;
			double l0 = trajectory_len - n_heads * arrow_head_step;
			double current_len = 0;
			size_t k=0;
			std::vector<std::array<double,4>> heads(n_heads);

			/* If the trajectory is long enough, calculate the envelope */
			const size_t NT = 2*len-2;
			std::vector<std::pair<double,double>> polygon(NT);
			polygon[0] = trajectory[t.j0];
			polygon[NT/2] = trajectory[t.j1-1];
			for (size_t j=t.j0+1; j<t.j1-1; ++j){
				/* Obtain width at current position: */
				double w = widths[j];

				/* Approximate symmetric direction at current index, rotated
				 * by 90°: */
				double dx2 = trajectory[j+1].second - trajectory[j-1].second;
				double dy2 = -(trajectory[j+1].first - trajectory[j-1].first);
				double invnorm = w/std::sqrt(dx2*dx2+dy2*dy2);
				dx2 *= invnorm;
				dy2 *= invnorm;

				/* Obtain points: */
				double x0 = trajectory[j].first + dx2;
				double y0 = trajectory[j].second + dy2;
				double x1 = trajectory[j].first - dx2;
				double y1 = trajectory[j].second - dy2;

				/* Add to point set: */
				points.add(trajectory[j], i+1);

				/* Save to polygon: */
				size_t p = j - t.j0; // Index in polygon
				polygon[p].first = x0;
				polygon[p].second = y0;
				polygon[NT-p].first = x1;
				polygon[NT-p].second = y1;

				/* See if we want to add another arrow head: */
				double dx = trajectory[j].first - trajectory[j-1].first;
				double dy = trajectory[j].second - trajectory[j-1].second;
				double len = std::sqrt(dx*dx + dy*dy);
				invnorm = 1.0 / len;
				dx *= invnorm;
				dy *= invnorm;
				current_len += len;
				if (current_len >= l0 + k*arrow_head_step){
					if (k == n_heads){
						throw std::runtime_error("Index j==n_heads out of bounds!");
					}

					heads[k] = {trajectory[j].first,trajectory[j].second, dx, dy};
					++k;
				}
			}

			/* Now add the envelope to list: */
			trajectories.push_back(polygon);
			all_heads.insert(all_heads.end(), heads.begin(), heads.end());
		}
	}

	return {trajectories, all_heads};
}

} // End namespace
