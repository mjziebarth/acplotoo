

#include <grid.hpp>
#include <field.hpp>
#include <utility>

namespace acplotoo {

typedef Field<std::pair<double,double>> Vectorfield;

typedef Field<double> Scalarfield;

/* 65536 trajectories seems like a reasonable safe upper bound to
 * what could possibly be visually pleasing.
 * This also corresponds to a starting grid of 256x256 points. */
typedef unsigned short tid_t;

typedef Field<tid_t> Mask;

std::pair<std::list<std::vector<std::pair<double,double>>>,
          std::vector<std::array<double,4>>>
    streamplot_polygons(const std::vector<std::pair<double,double>>& start,
                        const Vectorfield& velocity, const Grid& velocity_grid,
                        const Scalarfield& weights,
                        const std::pair<size_t,size_t>& mask_size,
                        double min_len, double max_len, double step_len_min,
                        double arrow_head_step, double collision_radius,
                        size_t max_steps, bool forward,
                        bool backward, double tol, size_t tile_history_size);

}
