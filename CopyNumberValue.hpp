#ifndef genomic_CopyNumberValue_h
#define genomic_CopyNumberValue_h

#include <cmath>

#include "typedefs.h"

template <typename a_type, typename b_type>
struct allele_specific
{
	a_type a;
	b_type b;
};

// sum of the absolute differences in each element
template <typename allele_specific_type> inline
allele_specific_type operator-(const allele_specific_type& x, const allele_specific_type& y) {
	return std::abs(x.a - y.a) + std::abs(x.b - y.b);
}

template <typename allele_specific_type> inline
bool operator!=(const allele_specific_type& x, const allele_specific_type& y) {
	return !( eq(x.a, y.a) && eq(y.b, y.b) );
}

template <typename allele_specific_type> inline
bool operator==(const allele_specific_type& x, const allele_specific_type& y) {
	return eq(x.a, y.a) && eq(x.b, y.b);
}

template <typename allele_specific_type> inline
std::istream& operator>>(std::istream& is, allele_specific_type& x)  {
	return is >> x.a >> ' ' >> x.b;
}

template <typename allele_specific_type>
inline bool eq(const allele_specifc_type& x, const allele_specific_type& y) {
	return eq(x.a, y.a) && eq(x.b, y.b);
}


struct AlleleSpecificCopyNumberValue
{
	CopyNumberValue a, b;
};

struct AlleleSpecificIntegerCopyNumberValue
{
	IntegerCopyNumberValue a, b;
};

inline CopyNumberValue operator-(const AlleleSpecificCopyNumberValue& x, const AlleleSpecificCopyNumberValue& y) {
	return std::abs(x.a - y.a) + std::abs(x.b - y.b);
}
inline bool operator!=(const AlleleSpecificCopyNumberValue& x, const AlleleSpecificCopyNumberValue& y) {
	return !(eq(x.a, y.a) && eq(x.b, y.b));
}
inline bool operator==(const AlleleSpecificCopyNumberValue& x, const AlleleSpecificCopyNumberValue& y) {
	return eq(x.a, y.a) && eq(x.b, y.b);
}
inline std::istream& operator>>(std::istream& is, AlleleSpecificCopyNumberValue& x)  {
	return is >> x.a >> ' ' >> x.b;
}

inline IntegerCopyNumberValue operator-(const AlleleSpecificIntegerCopyNumberValue& x, const AlleleSpecificIntegerCopyNumberValue& y) {
	return std::abs(x.a - y.a) + std::abs(x.b - y.b);
}
inline bool operator!=(const AlleleSpecificIntegerCopyNumberValue& x, const AlleleSpecificIntegerCopyNumberValue& y) {
	return !(x.a == y.a && y.b == y.b);
}
inline bool operator==(const AlleleSpecificIntegerCopyNumberValue& x, const AlleleSpecificIntegerCopyNumberValue& y) {
	return x.a == y.a && x.b == y.b;
}
inline std::istream& operator>>(std::istream& is, AlleleSpecificIntegerCopyNumberValue& x)  {
	return is >> x.a >> ' ' >> x.b;
}

inline bool eq(const AlleleSpecificIntegerCopyNumberValue& a, const AlleleSpecificIntegerCopyNumberValue& b) {
	return a == b;
}

inline bool eq(const AlleleSpecificCopyNumberValue& a, AlleleSpecificCopyNumberValue& b) {
	return a == b;
}

#endif