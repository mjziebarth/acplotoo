

#include <grid.hpp>
#include <field.hpp>
#include <utility>

namespace acplotoo {

typedef Field<std::pair<double,double>> Vectorfield;

typedef Field<size_t> Mask;

std::list<std::vector<std::pair<double,double>>>
    streamplot_polygons(const std::vector<std::pair<double,double>>& start,
                        const Vectorfield& velocity, const Grid& velocity_grid,
                        const std::pair<size_t,size_t>& mask_size,
                        double min_len, double max_len, size_t max_steps,
                        bool forward, bool backward, double tol);

}
