/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#ifndef compression_settings_h
#define compression_settings_h

#include <iostream>
#include "params.h"
#include "defs.h"

class CompSettings 
{
public:
    uint64_t vec_len;
    uint32_t max_no_vec_in_block;
    uint32_t ploidy;
    uint32_t max_depth;
    uint32_t ones_ranges;
    uint32_t n_samples;
    
    uint32_t n_vec_history_parts;
    uint32_t n_vec_history_vecs;
    
    uint32_t bit_size_id_match_pos_diff;
    uint32_t bit_size_id_copy_pos_diff;
    uint32_t bit_size_id;
    uint32_t bit_size_match_len;
    uint32_t bit_size_run_len = 0;
    uint32_t bit_size_literal;
    uint32_t bit_size_ones_goup;

    char bits_used(unsigned int n);
    
    CompSettings()
    {
        n_samples = 0;
        vec_len =  0;
        max_no_vec_in_block = 0;
        ploidy = 0;
        
        max_depth = 0;
        ones_ranges = 0;
        
        bit_size_id_match_pos_diff = 0;
        bit_size_id_copy_pos_diff  = 0;
        n_vec_history_parts = 0;
        n_vec_history_vecs = 0;
        
        bit_size_id = 0;
        bit_size_match_len = 0;
        bit_size_run_len = 0;
        bit_size_literal = 8;
        bit_size_ones_goup = 0;
    }
    
    CompSettings(Params params, uint32_t _n_samples)
    {
        n_samples = _n_samples;
        vec_len =  (n_samples *  params.ploidy) / 8 + (((n_samples * params.ploidy) % 8)?1:0) ;
        max_no_vec_in_block = params.var_in_block*2;
        ploidy = params.ploidy;
        
        max_depth = params.max_depth;
        ones_ranges = params.ones_ranges;
        
        bit_size_id_match_pos_diff = params.max_bit_size_id_match_pos_diff;
        bit_size_id_copy_pos_diff  = params.max_bit_size_id_copy_pos_diff;
        n_vec_history_parts = (1 << bit_size_id_match_pos_diff) ;
        n_vec_history_vecs = (1 << bit_size_id_copy_pos_diff);
        
        bit_size_id = (uint32_t)log2(max_no_vec_in_block) + 1;
        bit_size_match_len = (uint32_t)log2(vec_len) + 1;
        bit_size_run_len = (uint32_t)log2(vec_len) + 1;
        bit_size_literal = 8;
        bit_size_ones_goup = bits_used(ones_ranges - 1);
    }
};

#endif /* compression_settings_h */
