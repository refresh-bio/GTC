/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#ifndef _PARAMS_H
#define _PARAMS_H

#include <string>
#include <stdlib.h>
#include "defs.h"

struct Params{
    task_type task;
    file_type in_type, out_type;
    uint32_t max_depth, ones_ranges, max_MB_memory, max_bit_size_id_match_pos_diff, max_bit_size_id_copy_pos_diff;
    std::string in_file_name, in_ind_file;
    std::string arch_name;
    uint32_t ploidy;
    
    std::string samples;
    std::string range;
    
    std::string out_name;
    
    uint32_t sample_to_dec;
    uint32_t var_to_dec;
    uint32_t n_threads;
    uint32_t var_in_block;
    uint32_t records_to_process;
    
    char compression_level, mode;
    bool MB_memory;
    bool preprocessVCFonly;
    bool input_is_bit_vector;
    
    bool dec_single_sample, dec_single_var, out_permuted_bv;
    
    bool out_AC_AN, out_genotypes;
    
    uint32_t minAC, maxAC;
    float minAF, maxAF;
    
    Params()
    {
        in_type = VCF;
        out_type = VCF;
        max_depth = 100;
        ploidy = 2;
        n_threads = 2;
        var_in_block = PART_SIZE;
        arch_name = "archive";
        out_name = "";
        range = "";
        samples = "";
        ones_ranges = 8;
        MB_memory = true;  // Remember some of the decoded vectors
        max_MB_memory = 0;  // 0 means no limit (if MB_memory == true)
        compression_level = '1';
        max_bit_size_id_match_pos_diff = 8; //(1 << 17) - 1;
        max_bit_size_id_copy_pos_diff  = 17; //(1 << 17) - 1;
        
        preprocessVCFonly = false;
        input_is_bit_vector = false;
        
        dec_single_sample = false;
        dec_single_var = false;
        out_permuted_bv = false;
        
        out_AC_AN = false;
        out_genotypes = true;
        records_to_process = UINT32_MAX;
        mode = '\0';
        
        minAC = 0;
        maxAC = INT32_MAX;
        
        minAF = 0;
        maxAF = 1;
    }
};

#endif




