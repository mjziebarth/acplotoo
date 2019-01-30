#include <pointset.hpp>

#include <exception>

namespace acplotoo {

PointSet::PointSet()
{}

void PointSet::add(std::pair<double,double> point, unsigned short id)
{
	x_map.emplace(point.first, value_t(point.second, id));
}

std::vector<std::tuple<double,double,unsigned short>>
PointSet::query_in_range(double x, double y, double rmax) const
{
	/* Sanity check on radius: */
	if (rmax < 0.0){
		throw std::domain_error("Radius must be non-negative.");
	}
	const double rmax2 = rmax*rmax;

	/* Find data in x range: */
	auto lower = x_map.lower_bound(x-rmax);
	auto upper = x_map.upper_bound(x+rmax);

	/* Iterate over range and for each element check distance: */
	std::vector<std::tuple<double,double,unsigned short>> result;
	for (auto it=lower; it != upper; ++it){
		double dx = x - it->first;
		double dy = y - it->second.first;
		if (dx*dx + dy*dy <= rmax2){
			/* Element within range. */
			result.emplace_back(it->first, it->second.first, it->second.second);
		}
	}

	/* Return vector: */
	return result;
}

void PointSet::remove_all(unsigned short id)
{
	/* Removes all entries of id. */
	for (auto it = x_map.begin(); it != x_map.end();){
		if (it->second.second == id){
			it = x_map.erase(it);
		} else {
			++it;
		}
	}
}

} // END NAMESPACE
