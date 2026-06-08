#ifndef cna_global_h
#define cna_global_h

#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <limits>
#include <cmath>

#include "logging.hpp"

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>

#include "typedefs.h"

#include "config.h"

namespace cna {

template <typename T>
inline typename boost::enable_if<boost::is_arithmetic<T>, T>::type absdiff(T a, T b) {
	return (a > b) ? a - b : b - a;
}

template <typename T>
inline typename boost::enable_if<boost::is_arithmetic<T>, bool>::type eq(T a, T b, T epsilon=std::numeric_limits<T>::epsilon()) {
	return absdiff(a, b) <= epsilon;
}

template <typename T>
inline typename boost::enable_if<boost::is_arithmetic<T>, bool>::type neq(T a, T b, T epsilon=std::numeric_limits<T>::epsilon()) {
	return absdiff(a, b) > epsilon;
}

// Number of human autosomes
extern chromid nAutosomes;
// Number of human chromosomes: 22 autosomes + 2 sex chromosomes
extern chromid nChromosomes;

namespace data
{
	enum Type {
		invalid, generic, raw, segmented, raw_ref, segmented_ref, raw_ascn, segmented_ascn, raw_lrrbaf
	};
};

namespace mapping
{
	/* singleton */
	class ChromosomesMap
	{
	private:
		typedef std::map<std::string, chromid> chr2index;
		typedef std::map<chromid, std::string> index2chr;
		chr2index index;
		index2chr chr;
	public:
		ChromosomesMap() {
			
			// Map numeric values of chromosomes and alternatives prefixed with chr
			char s[6];
			for (chromid i = 1; i <= cna::nChromosomes; ++i) {
				sprintf(s, "%d", i); index[s] = i;
				index[s] = i;
			}
			
			// Generate reverse mapping (while the mapping is still one-to-one)
			chr2index::const_iterator it, end = index.end();
			for (it = index.begin(); it != end; ++it) {
				chr.insert(std::make_pair(it->second, it->first));
			}
			
			// Map special cases
			index["X"] = 23;
			index["Y"] = 24;
			chr[23] = "X";
			chr[24] = "Y";
			
			// Map alternative names for chromosomes prefixed with chr
			for (chromid i = 1; i <= cna::nChromosomes; ++i) {
				sprintf(s, "chr%d", i);
				index[s] = i;
			}
			index["chrX"] = 23;
			index["chrY"] = 24;
		}	
		
		chromid operator[] (const std::string& chr) {
			return index[chr];
			/*
			chr2index::const_iterator it = index.find(chr);
			if (it == index.end()) {
				return cna::nChromosomes+1;
			}
			return it->second;
			*/
		}
		
		std::string operator[] (chromid index) {
			return chr[index];
			/*
			index2chr::const_iterator it = chr.find(index);
			if (it == chr.end()) {
				return "";
			}
			return it->second;
			*/
		}
		
		void print() {
			chr2index::const_iterator end = index.end();
			for (std::map<std::string, chromid>::iterator it = index.begin(); it != end; ++it) {
				std::cout << it->first << " -> " << it->second << std::endl;
			}
		}
	};
	extern ChromosomesMap chromosome;
	
	/* singleton */
	class ExtensionMap
	{
	private:
		typedef std::map<std::string, data::Type> ext2type;
		typedef std::map<data::Type, std::string> type2ext;
		ext2type type;
		type2ext ext;
	public:
		ExtensionMap() {
			type["cn"] = data::raw;
			type["seg"] = data::segmented;
			type["cnref"] = data::raw_ref;
			type["segref"] = data::segmented_ref;
			type["cnas"] = data::raw_ascn;
			type["segas"] = data::segmented_ascn;
			type["lrrbaf"] = data::raw_lrrbaf;
			
			// create reverse mapping
			ext2type::const_iterator it, end = type.end();
			for (it = type.begin(); it != end; ++it) {
				ext.insert(std::make_pair(it->second, it->first));
			}
		}
		data::Type operator[] (const std::string& ext) {
			ext2type::const_iterator it = type.find(ext);
			if (it == type.end()) {
				return data::invalid;
			} else {
				return it->second;
			}
		}
		const std::string& operator[] (data::Type type) {
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
	bool pair(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) {
		return a.first < b.first;
	}
}

namespace name
{
	std::string common(const std::string& a, const std::string& b);
	std::string fileext(const std::string& s);
	std::string filestem(const std::string& s);
	std::string filepath(const std::string& s);
}

} // namespace cna

#endif
