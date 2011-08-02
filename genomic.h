#ifndef _genomic_h_
#define _genomic_h_

#if genomic_DEBUG == 1
#define trace(...) printf(__VA_ARGS__)
#else
#define trace(...)
#endif

#include "config.h"
#include "global.h"
#include "Sample.h"

#endif