#ifndef genomic_CopyNumberValue_h
#define genomic_CopyNumberValue_h

#include <cmath>
#include <limits>

#include "typedefs.h"

struct AlleleSpecificCopyNumberValue
{
	CopyNumberValue a, b;
};

struct AlleleSpecificIntegerCopyNumberValue
{
	IntegerCopyNumberValue a, b;
};

template <typename T>
inline bool eq(T a, T b, T epsilon=std::numeric_limits<T>::epsilon()) {
	// essentially equal
	//return abs(a - b) <= ( (abs(a) > abs(b) ? abs(b) : abs(a)) * epsilon );
	return std::abs(a - b) <= epsilon;
}

template <typename T>
inline bool neq(T a, T b, T epsilon=std::numeric_limits<T>::epsilon()) {
	// essentially equal
	//return abs(a - b) > ( (abs(a) > abs(b) ? abs(b) : abs(a)) * epsilon );
	return std::abs(a - b) > epsilon;
}

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
	return is >> x.a >> x.b;
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
	return is >> x.a >> x.b;
}

inline bool eq(const AlleleSpecificIntegerCopyNumberValue& a, const AlleleSpecificIntegerCopyNumberValue& b) {
	return a == b;
}

inline bool eq(const AlleleSpecificCopyNumberValue& a, AlleleSpecificCopyNumberValue& b) {
	return a == b;
}

#endif