/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#ifndef my_vcf_h
#define my_vcf_h

#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "htslib/kstring.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/hfile.h"
#include "htslib/khash_str2int.h"


#include "htslib/khash.h"
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

#include "htslib/kseq.h"
KSTREAM_DECLARE(gzFile, gzread)


#define hts_expand00(type_t, n, m, ptr) if ((n) > (m)) { \
int t = (m); (m) = (n); kroundup32(m); \
(ptr) = (type_t*)realloc((ptr), (m) * sizeof(type_t)); \
memset(((type_t*)ptr)+t,0,sizeof(type_t)*((m)-t)); \
}


int bcf_update_genotypes_fast(kstring_t &s, bcf1_t *line);

#endif /* my_vcf_h */
