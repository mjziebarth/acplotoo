/* field.hpp */

#ifndef ACPLOTOO__FIELD_HPP
#define ACPLOTOO__FIELD_HPP

#include <vector>
#include <utility>
#include <grid.hpp>
#include <exception>
#include <algorithm>

namespace acplotoo {

/* Class declaration: */
template<typename T>
class Field {
	public:
		typedef Grid::index_t index_t;
		/* M: Size of y coordinates.*/
		Field(size_t M, size_t N, T default_value=T(), bool y_first=true);

		/* Accessing elements (performs bounds checking): */
		T& operator[](const index_t& index);
		const T& operator[](const index_t& index) const;

		/* For easy Cython interfacing: */
		void set(size_t i, size_t j, const T& val);

		size_t nx() const;
		size_t ny() const;

		/* Replace all elements of type: */
		void replace_all(const T& which, const T& by);

	private:
		size_t M,N;
		size_t y_first;

		std::vector<std::vector<T>> data;

};


/* Template code: */



template<typename T>
Field<T>::Field(size_t M, size_t N, T default_value, bool y_first) : M(M), N(N), y_first(y_first),
    data((y_first) ? N : M)
{
	/* Initialize matrix as a vector of vectors: */
	if (y_first){
		for (size_t i=0; i<N; ++i){
			data[i].resize(M, default_value);
		}
	} else {
		for (size_t i=0; i<M; ++i){
			data[i].resize(N, default_value);
		}
	}
}


template<typename T>
T& Field<T>::operator[](const index_t& index)
{
	const size_t i = index.first;
	const size_t j = index.second;

	/* Be safe! */
	if (i > M || j > N)
		throw std::out_of_range("Field: Tried to access element out of range!");

	if (y_first)
		return data[j][i];

	return data[i][j];
}

template<typename T>
const T& Field<T>::operator[](const index_t& index) const
{
	const size_t i = index.first;
	const size_t j = index.second;

	/* Be safe! */
	if (i > M || j > N)
		throw std::out_of_range("Field: Tried to access element out of range!");

	if (y_first)
		return data[j][i];

	return data[i][j];
}

template<typename T>
size_t Field<T>::ny() const
{
	return N;
}

template<typename T>
size_t Field<T>::nx() const
{
	return M;
}

template<typename T>
void Field<T>::replace_all(const T& which, const T& by)
{
	for (std::vector<T>& v : data){
		std::replace(v.begin(), v.end(), which, by);
	}
}





} // End namespace

#endif
