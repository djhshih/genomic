#ifndef genomic_Chromosome_h
#define genomic_Chromosome_h

#include <vector>

#include <boost/type_traits/is_pointer.hpp>
#include <boost/type_traits/remove_pointer.hpp>

//#include "Tree.hpp"


template <typename V> class Segment;

template <typename T>
class Chromosome
{
public:
	typedef T DataType;
private:
	
public:
	size_t index;
	Chromosome(chromid chromIndex) : index(chromIndex) {}
};

template <typename T>
class LinearChromosome : Chromosome<T>
{
public:
	typedef T DataType;
	typedef typename std::vector<T>::iterator iterator;
	typedef typename std::vector<T>::const_iterator const_iterator;
private:
	std::vector<T> items;
public:
	LinearChromosome(chromid chromIndex) : Chromosome<T>(chromIndex) {}
	LinearChromosome(const LinearChromosome& chr) 
	: Chromosome<T>(chr.index), items(chr.items) {
		// clone copy if T is a pointer type
		duplicate( typename boost::is_pointer<T>::type() );
	}
	
	LinearChromosome<T>& operator=(const LinearChromosome<T>& chr) {
		LinearChromosome tmp(chr);
		swap(tmp);
		return *this;
	}
	LinearChromosome<T>& swap(LinearChromosome<T>& chr) {
		items.swap(chr.items);
	}
	~LinearChromosome() {
	}
	T& at(size_t i) {
		return items[i];
	}
	T& operator[](size_t i) {
		return items[i];
	}
	const size_t size() const {
		return items.size();
	}
	iterator begin() {
		return items.begin();
	}
	iterator end() {
		return items.end();
	}
	const_iterator begin() const {
		return items.begin();
	}
	const_iterator end() const {
		return items.end();
	}
	void push_back(const T& item) {
		items.push_back(item);
	}
	void pop_back() {
		items.pop_back();
	}
	void clear() {
		// deallocate memory if T is a pointe type
		free( typename boost::is_pointer<T>::type() );
		items.clear();
	}
	void resize(size_t n) {
		items.resize(n);
	}
	
private:
	
	void duplicate(const typename boost::true_type& t) {
		for (chromid i = 0; i < items.size(); ++i) {
			items[i] = new typename boost::remove_pointer<T>::type(*(items[i]));
		}
	}
	void duplicate(const typename boost::false_type& t) {}
	
	void free(const typename boost::true_type&) {
		for (chromid i = 0; i < items.size(); ++i) {
			delete items[i];
		}
	}
	void free(const typename boost::false_type&) {}
	
};


/*
template <typename T>
class KDTreeChromosome : Chromosome<T>
{
private:
	tree::KDTree< tree::Container< Segment<T> > > items;
public:
	KDTreeChromosome(chromid chromIndex) {
	}
};
*/

#endif