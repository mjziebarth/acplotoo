

#ifndef POINTSET_HPP
#define POINTSET_HPP

#include <map>
#include <vector>
#include <utility>
#include <tuple>


namespace acplotoo {

class PointSet {
	public:
		PointSet();

		/* Add an ID at a point: */
		void add(std::pair<double,double> point, unsigned short id);

		/* Query all IDs and the associated coordinates
		 * in range of a point: */
		std::vector<std::tuple<double,double,unsigned short>>
		    query_in_range(double x, double y, double rmax) const;

		/* Remove all by id: */
		void remove_all(unsigned short id);

	private:
		typedef std::pair<double, unsigned short> value_t;
	
		std::multimap<double, value_t> x_map;
};

} // end namespace

#endif
