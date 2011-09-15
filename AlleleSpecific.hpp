#ifndef genomic_AlleleSpecific_h
#define genomic_AlleleSpecific_h

#include <cmath>
#include <limits>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_base_of.hpp>

#include "typedefs.h"


#define ENABLE_IF_ALLELE_SPECIFIC typename boost::enable_if<boost::is_base_and_derived<allele_specific_base, allele_specific_type>, allele_specific_type>

struct allele_specific_base {};

template <typename a_type, typename b_type>
struct allele_specific : allele_specific_base
{
	a_type a;
	b_type b;
};

// sum of the absolute differences in each element
template <typename allele_specific_type> inline
allele_specific_type operator-(const allele_specific_type& x, const ENABLE_IF_ALLELE_SPECIFIC::type& y) {
	return std::abs(x.a - y.a) + std::abs(x.b - y.b);
}

template <typename allele_specific_type> inline
bool operator<=(const allele_specific_type& x, float z) {
	return x.a <= z;
}

template <typename allele_specific_type> inline
bool operator>=(const allele_specific_type& x, float z) {
	return x.a >= z;
}

inline
bool operator<=(const alleles_rcn& x, float z) {
	return (x.a + x.b) <= z;
}

inline
bool operator>=(const alleles_rcn& x, float z) {
	return (x.a + x.b) >= z;
}

template <typename allele_specific_type> inline
bool operator!=(const allele_specific_type& x, const ENABLE_IF_ALLELE_SPECIFIC::type& y) {
	return !( eq(x.a, y.a) && eq(y.b, y.b) );
}

template <typename allele_specific_type> inline
bool operator==(const allele_specific_type& x, const ENABLE_IF_ALLELE_SPECIFIC::type& y) {
	return eq(x.a, y.a) && eq(x.b, y.b);
}

template <typename allele_specific_type> inline
std::istream& operator>>(std::istream& is, ENABLE_IF_ALLELE_SPECIFIC::type& x)  {
	return is >> x.a >> ' ' >> x.b;
}

template <typename allele_specific_type>
inline bool eq(const allele_specific_type& x, const ENABLE_IF_ALLELE_SPECIFIC::type& y) {
	return eq(x.a, y.a) && eq(x.b, y.b);
}

template <typename allele_specific_type>
inline bool neq(const allele_specific_type& x, const ENABLE_IF_ALLELE_SPECIFIC::type& y) {
	return neq(x.a, y.a) || neq(x.b, y.b);
}

#endif