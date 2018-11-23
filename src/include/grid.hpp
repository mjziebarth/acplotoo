/* grid.hpp */

#ifndef ACPLOTOO__GRID_HPP
#define ACPLOTOO__GRID_HPP

#include <utility>
#include <array>
#include <string>

namespace acplotoo {

/* Point class */
class Point {
	public:
		Point(double x, double y);

		double x() const;
		double y() const;

		operator std::pair<double,double>() const;

		std::string to_string() const;

	private:
		double x_, y_;
};


/* Rect class (template) */
enum vertex_t {
	BOT_LEFT, BOT_RIGHT, TOP_LEFT, TOP_RIGHT
};

template<typename T>
class Rect {
	public:
		Rect(T bot_left, T bot_right, T top_left, T top_right);

		T operator[](vertex_t) const;

	private:
		T bot_left, bot_right, top_left, top_right;
};

template<typename T>
Rect<T>::Rect(T bot_left, T bot_right, T top_left, T top_right) : bot_left(bot_left),
   bot_right(bot_right), top_left(top_left), top_right(top_right)
{
}

template<typename T>
T Rect<T>::operator[](vertex_t v) const
{
	if (v == BOT_LEFT){
		return bot_left;
	} else if (v == BOT_RIGHT){
		return bot_right;
	} else if (v == TOP_LEFT){
		return top_left;
	}
	return top_right;
}




/* Grid class */
class Grid {
	public:
		typedef std::pair<size_t,size_t> index_t;

		Grid(double x0, double x1, double y0, double y1, size_t nx, size_t ny);

		Grid resample(size_t nx, size_t ny) const;

		/* Access grid element at (i,j): */
		Point operator[](const index_t& index) const;

		/* Access closest grid point: */
		index_t closest(double x, double y) const;

		/* Check if element is contained: */
		bool contains(double x, double y) const;

		/* Obtain indices of the nodes that make up the cell this
		 * point is in: */
		Rect<index_t> cell(double x, double y) const;

		/* Check if the space covered is the same as that of another
		 * grid: */
		bool same_rect(const Grid& other) const;

		/* Return size: */
		std::pair<size_t,size_t> shape() const;

		/* Compare grid: */
		bool operator==(const Grid& other) const;

	private:
		double x0, x1, width, dx, y0, y1, height, dy;
		size_t nx, ny;
};

} // Namespace

#endif
