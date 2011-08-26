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

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>

#include "typedefs.h"

#include "config.h"
#if genomic_DEBUG == 1
#define trace(...) printf(__VA_ARGS__)
#else
#define trace(...)
#endif


#define ENABLE_IF_ARITHMETIC typename boost::enable_if<boost::is_arithmetic<T>, T>

template <typename T>
inline bool eq(T a, T b, ENABLE_IF_ARITHMETIC::type epsilon=std::numeric_limits<T>::epsilon()) {
	// essentially equal
	//return abs(a - b) <= ( (abs(a) > abs(b) ? abs(b) : abs(a)) * epsilon );
	return std::abs(a - b) <= epsilon;
}

template <typename T>
inline bool neq(T a, T b, ENABLE_IF_ARITHMETIC::type epsilon=std::numeric_limits<T>::epsilon()) {
	// essentially equal
	//return abs(a - b) > ( (abs(a) > abs(b) ? abs(b) : abs(a)) * epsilon );
	return std::abs(a - b) > epsilon;
}


// Number of human autosomes
extern chromid nAutosomes;
// Number of human chromosomes: 22 autosomes + 2 sex chromosomes
extern chromid nChromosomes;

namespace data
{
	enum Type {
		invalid, generic, raw, segmented, raw_ascn, segmented_ascn, dchip, cnag, picnic, penncnv
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
	extern ChromosomesMap chromosome;
	
	/* singleton */
	class ExtensionMap
	{
	private:
		typedef map<string, data::Type> ext2type;
		typedef map<data::Type, string> type2ext;
		ext2type type;
		type2ext ext;
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
			
			// create reverse mapping
			ext2type::const_iterator it, end = type.end();
			for (it = type.begin(); it != end; ++it) {
				ext.insert( make_pair(it->second, it->first) );
			}
		}
		data::Type operator[] (const string& ext) {
			ext2type::const_iterator it = type.find(ext);
			if (it == type.end()) {
				return data::invalid;
			} else {
				return it->second;
			}
		}
		const string& operator[] (data::Type type) {
			type2ext::const_iterator it = ext.find(type);
			if (it == ext.end()) {
				return ext[data::invalid];
			} else {
				return it->second;
			}
		}
	};
	extern ExtensionMap extension;
	
}

namespace compare
{
	template <typename T1, typename T2>
	bool pair(const pair<T1, T2>& a, const pair<T1, T2>& b) {
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