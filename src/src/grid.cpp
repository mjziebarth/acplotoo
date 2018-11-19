
#include <grid.hpp>
#include <exception>
#include <cmath>

namespace acplotoo {


/* Point: */
Point::Point(double x, double y) : x_(x), y_(y)
{
}

Point::operator std::pair<double,double>() const
{
	return std::pair<double,double>(x_,y_);
}

double Point::x() const
{
	return x_;
}

double Point::y() const
{
	return y_;
}



/* Grid: */

Grid::Grid(double x0, double x1, double y0, double y1, size_t nx, size_t ny)
	: x0(x0), y0(y0), x1(x1), y1(y1), width(x1-x0), height(y1-y0), nx(nx), ny(ny),
	  dx(width / nx), dy(height / ny)
{
	/* Ensure limits are correct: */
	if (x0 >= x1)
		throw std::runtime_error("Wrong order of x limits or empty grid!");
	if (y0 >= y1)
		throw std::runtime_error("Wrong order of y limits or empty grid!");
	if (nx == 0 || ny == 0)
		throw std::runtime_error("Empty grid!");
}

Grid Grid::resample(size_t nx, size_t ny) const
{
	return Grid(x0, x1, y0, y1, nx, ny);
}


Point Grid::operator[](const Grid::index_t& index) const
{
	const size_t i=index.first;
	const size_t j=index.second;
	/* Range check: */
	if (i >= nx || j >= ny)
		throw std::out_of_range("Grid: Tried to access element out of range!");
	if (i == nx-1){
		if (j == ny-1){
			return Point(x1,y1);
		}
		return Point(x1, y0 + j*dy);
	} else if (j == ny-1) {
		return Point(x0 + i*dx, y1);
	}
	return Point(x0 + i*dx, y0 + j*dy);
}


Grid::index_t Grid::closest(double x, double y) const
{
	size_t i,j;
	/* x logic */
	if (x <= x0)
		i = 0;
	else if (x >= x1)
		i = nx-1;
	else
		i = std::lround(nx * (x-x0) / width);
	/* y logic: */
	if (y <= y0)
		j = 0;
	else if (y >= y1)
		j = ny-1;
	else
		j = std::lround(ny * (y-y0) / height);

	return index_t(i,j);
}


bool Grid::contains(double x, double y) const
{
	return (x >= x0) && (x <= x1) && (y >= y0) && (y <= y1);
}

Rect<Grid::index_t> Grid::cell(double x, double y) const
{
	/* Make sure that point is contained: */
	if (!contains(x,y))
		throw std::runtime_error("Grid::cell(): Point not contained!");


	/* Obtain the rect that is described by respective flooring and
	 * ceiling of x and y coordinate. Use min/max to make sure that
	 * index is in bounds. */
	double sx = (x-x0) / dx;
	double sy = (y-y0) / dy;
	size_t i0 = std::min(std::max<size_t>(std::floor(sx), 0), nx-1);
	size_t i1 = std::min(std::max<size_t>(std::ceil(sx), 0), nx-1);
	size_t j0 = std::min(std::max<size_t>(std::floor(sy), 0), ny-1);
	size_t j1 = std::min(std::max<size_t>(std::ceil(sy), 0), ny-1);

	/* Create rect: Rect(bot_left, bot_right, top_left, top_right)*/
	typedef Rect<Grid::index_t> rect_t;
	return Rect<index_t>(index_t(i0,j0), index_t(i1,j0),
	                     index_t(i0,j1), index_t(i1,j1));
}

bool Grid::same_rect(const Grid& other) const
{
	return (x0 == other.x0) && (x1 == other.x1) && (y0 == other.y0)
	        && (y1 == other.y1);
}

std::pair<size_t,size_t> Grid::shape() const
{
	return std::pair<size_t,size_t>(nx,ny);
}


} // NAMESPACE
