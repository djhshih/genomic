#ifndef genomic_typedefs_h
#define genomic_typedefs_h

// chromosome index
typedef unsigned int chromid;

// genomic position
typedef unsigned long position;

// difference in genomic position
typedef signed long position_diff;

// copy number value
typedef unsigned int cnvalue;

// relative copy number value
typedef int rcnvalue;

// raw total signal
typedef float rvalue;

// log R ratio
typedef float lrratio;

// B allele frequency
typedef float bafreq;

// allele-specific type template
template <typename a_type, typename b_type> struct allele_specific;

// allele specific number copy information, available as LRR and BAF
typedef allele_specific<lrratio, bafreq> alleles_raw;

// total copy number value and minor copy number value
typedef allele_specific<cnvalue, cnvalue> alleles_cn;

// relative copy number value for either allele
typedef allele_specific<rcnvalue, rcnvalue> alleles_rcn;


#endif
