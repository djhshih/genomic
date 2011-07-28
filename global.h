#ifndef genomic_global_h
#define genomic_global_h

#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <limits>

using namespace std;


typedef int IntegerCopyNumberValue;
typedef float CopyNumberValue;
typedef unsigned int chromid;

template <typename T>
inline bool eq(T a, T b, T epsilon=numeric_limits<T>::epsilon()) {
	// essentially equal
	//return abs(a - b) <= ( (abs(a) > abs(b) ? abs(b) : abs(a)) * epsilon );
	return abs(a - b) <= epsilon;
}

struct AlleleSpecificCopyNumberValue
{
	CopyNumberValue a, b;
};

inline CopyNumberValue operator-(const AlleleSpecificCopyNumberValue& x, const AlleleSpecificCopyNumberValue& y) {
	return abs(x.a - y.a) + abs(x.b - y.b);
}
inline bool operator!=(const AlleleSpecificCopyNumberValue& x, const AlleleSpecificCopyNumberValue& y) {
	return !(eq(x.a, y.a) && eq(x.b, y.b));
}
inline bool operator==(const AlleleSpecificCopyNumberValue& x, const AlleleSpecificCopyNumberValue& y) {
	return eq(x.a, y.a) && eq(x.b, y.b);
}
inline istream& operator>>(istream& is, AlleleSpecificCopyNumberValue& x)  {
	return is >> x.a >> x.b;
}

struct AlleleSpecificIntegerCopyNumberValue
{
	IntegerCopyNumberValue a, b;
};

inline IntegerCopyNumberValue operator-(const AlleleSpecificIntegerCopyNumberValue& x, const AlleleSpecificIntegerCopyNumberValue& y) {
	return abs(x.a - y.a) + abs(x.b - y.b);
}
inline bool operator!=(const AlleleSpecificIntegerCopyNumberValue& x, const AlleleSpecificIntegerCopyNumberValue& y) {
	return !(x.a == y.a && y.b == y.b);
}
inline bool operator==(const AlleleSpecificIntegerCopyNumberValue& x, const AlleleSpecificIntegerCopyNumberValue& y) {
	return x.a == y.a && x.b == y.b;
}
inline istream& operator>>(istream& is, AlleleSpecificIntegerCopyNumberValue& x)  {
	return is >> x.a >> x.b;
}

inline bool eq(const AlleleSpecificIntegerCopyNumberValue& a, const AlleleSpecificIntegerCopyNumberValue& b) {
	return a == b;
}

inline bool eq(const AlleleSpecificCopyNumberValue& a, AlleleSpecificCopyNumberValue& b) {
	return a == b;
}

typedef unsigned long position;
typedef signed long position_diff;

// Number of human autosomes
extern chromid nAutosomes;
// Number of human chromosomes: 22 autosomes + 2 sex chromosomes
extern chromid nChromosomes;

namespace data
{
	enum Type {
		generic, raw, segmented, raw_ascn, segmented_ascn, dchip, cnag, picnic, penncnv
	};
};

namespace mapping
{
	/* singleton */
	class ChromosomesMap
	{
	private:
		map<string, chromid> index;
	public:
		ChromosomesMap() {
			// Map numeric values of chromosomes and alternatives prefixed with chr
			char s[6];
			for (chromid i = 1; i <= nChromosomes; ++i) {
				sprintf(s, "%d", i); index[s] = i;
				sprintf(s, "chr%d", i);
				index[s] = i;
			}
			// Map alternative names for chromosomes
			index["X"] = 23;
			index["chrX"] = 23;
			index["Y"] = 24;
			index["chrY"] = 24;
		}	
		chromid operator[] (const string& chr) {
			return index[chr];
		}
		void print() {
			map<string, chromid>::const_iterator end = index.end();
			for (map<string, chromid>::iterator it = index.begin(); it != end; ++it) {
				cout << it->first << " -> " << it->second << endl;
			}
		}
	};
	static ChromosomesMap chromosome;
	
	/* singleton */
	class ExtensionMap
	{
	private:
		map<string, data::Type> type;
	public:
		ExtensionMap() {
			type["cn"] = data::raw;
			type["seg"] = data::segmented;
			type["cnas"] = data::raw_ascn;
			type["segas"] = data::segmented_ascn;
			type["dchip"] = data::dchip;
			type["cnag"] = data::cnag;
			type["picnic"] = data::picnic;
			type["penncnv"] = data::penncnv;
		}
		data::Type operator[] (const string& ext) {
			return type[ext];
		}
	};
	static ExtensionMap extension;
};

namespace compare
{
	template <typename T1, typename T2>
	bool pair(pair<T1, T2> a, pair<T1, T2> b) {
		return a.first < b.first;
	}
}

namespace name
{
	string common(const string& a, const string& b);
	string fileext(const string& s);
	string filestem(const string& s);
	string filepath(const string& s);
}

#endif