#ifndef _genomic_h_
#define _genomic_h_

#define DEBUG

#ifdef DEBUG
#define trace(...) printf(__VA_ARGS__)
#else
#define trace(...)
#endif

#include "global.h"
#include "Sample.h"

#endif