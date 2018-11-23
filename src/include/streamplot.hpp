

#include <grid.hpp>
#include <field.hpp>
#include <utility>

namespace acplotoo {

typedef Field<std::pair<double,double>> Vectorfield;

typedef Field<double> Scalarfield;

typedef Field<size_t> Mask;

std::pair<std::list<std::vector<std::pair<double,double>>>,
          std::vector<std::array<double,4>>>
    streamplot_polygons(const std::vector<std::pair<double,double>>& start,
                        const Vectorfield& velocity, const Grid& velocity_grid,
                        const Scalarfield& weights,
                        const std::pair<size_t,size_t>& mask_size,
                        double min_len, double max_len, double step_len_min,
                        double arrow_head_step, size_t max_steps, bool forward,
                        bool backward, double tol);

}
