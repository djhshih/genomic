#ifndef genomic_typedefs_h
#define genomic_typedefs_h

typedef unsigned int chromid;

typedef unsigned long position;
typedef signed long position_diff;

typedef int IntegerCopyNumberValue;
typedef float CopyNumberValue;

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


template <typename a_type, typename b_type> struct allele_specific;

// allele specific number copy information, available as LRR and BAF
typedef allele_specific<lrratio, bafreq> alleles_raw;

// total copy number value and minor copy number value
typedef allele_specific<cnvalue, cnvalue> alleles_cn;

// relative copy number value for either allele
typedef allele_specific<rcnvalue, rcnvalue> alleles_rcn;


#endif
