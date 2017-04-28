/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April

 It is a modification of vcf.c file by Heng Li (see note below)
 */

/*  vcf.c -- VCF/BCF API functions.

    Copyright (C) 2012, 2013 Broad Institute.
    Copyright (C) 2012-2015 Genome Research Ltd.
    Portions copyright (C) 2014 Intel Corporation.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include "my_vcf.h"



#define bit_array_size(n) ((n)/8+1)
#define bit_array_set(a,i)   ((a)[(i)/8] |=   1 << ((i)%8))
#define bit_array_clear(a,i) ((a)[(i)/8] &= ~(1 << ((i)%8)))
#define bit_array_test(a,i)  ((a)[(i)/8] &   (1 << ((i)%8)))

static inline uint8_t *bcf_unpack_fmt_core1(uint8_t *ptr, int n_sample, bcf_fmt_t *fmt)
{
    uint8_t *ptr_start = ptr;
    fmt->id = bcf_dec_typed_int1(ptr, &ptr);
    fmt->n = bcf_dec_size(ptr, &ptr, &fmt->type);
    fmt->size = fmt->n << bcf_type_shift[fmt->type];
    fmt->p = ptr;
    fmt->p_off  = ptr - ptr_start;
    fmt->p_free = 0;
    ptr += n_sample * fmt->size;
    fmt->p_len = ptr - fmt->p;
    return ptr;
}


// Based on function from htslib: bcf_update_format((hdr),(line),"GT",(gts),(n),BCF_HT_INT)
int bcf_update_genotypes_fast(kstring_t &str,  bcf1_t *line)
{
    bcf_fmt_t *fmt = NULL;

    line->n_fmt++;
    hts_expand00(bcf_fmt_t, line->n_fmt, line->d.m_fmt, line->d.fmt);
    fmt = &line->d.fmt[line->n_fmt-1];
    bcf_unpack_fmt_core1((uint8_t*)str.s, line->n_sample, fmt);
    line->d.indiv_dirty = 1;
    fmt->p_free = 0;

    line->unpacked |= BCF_UN_FMT;
    return 0;
}
