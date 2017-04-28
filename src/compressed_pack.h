/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#ifndef compressed_pack_h
#define compressed_pack_h

#include <stdio.h>
#include <iostream>
#include <string>
#include <sdsl/bit_vectors.hpp>

#include "defs.h"
#include "compression_settings.h"
#include "huffman.h"

#define MMAP

#ifdef MMAP
//#include <boost/interprocess/file_mapping.hpp>
//#include <boost/interprocess/mapped_region.hpp>
#include <cpp-mmf/memory_mapped_file.hpp>
#endif

class CompressedPack {
    friend class Decompressor;
    CompSettings s;
    uchar * buf = nullptr;
    
    // sdsl vectors and ranks
    sdsl::rrr_vector<> rrr_copy_bit_vector[2];  //for copy vectors (even and odd)
    sdsl::rrr_vector<> rrr_zeros_only_bit_vector[2];  //for zeros vectors (even and odd)

    sdsl::rrr_vector<>::rank_0_type rrr_rank_zeros_only_bit_vector_0;
    sdsl::rrr_vector<>::rank_1_type rrr_rank_zeros_only_bit_vector_1;
    sdsl::rrr_vector<>::rank_1_type rrr_rank_copy_bit_vector[2];
    
#ifdef MMAP
//    boost::interprocess::file_mapping *fm;
//    boost::interprocess::mapped_region *region;
    memory_mapped_file::read_only_mmf *fm;
#endif
    
    // huffmans
    CHuffman huf_match_diff_MSB;
    
    CHuffman * huf_literals = nullptr;
    CHuffman * huf_zeros_runs = nullptr;
    CHuffman * huf_ones_runs = nullptr;
    CHuffman * huf_match_lens = nullptr;
    CHuffman huf_flags;
    CHuffman huf_group_type;
    
    // metadata
    uint64_t no_vec;
    uint64_t no_copy;
    uint64_t no_non_copy;
    uint32_t used_bits_cp;
    uint32_t  used_bits_noncp;
    uint32_t  no_blocks;
    
    uint32_t  max_used_bits_litRunSize[MAX_NUMBER_OF_GROUP], used_bits_litRunSize[MAX_NUMBER_OF_GROUP] ;
    uint32_t minLitRunSize[MAX_NUMBER_OF_GROUP][MAX_LITERAL_RUN + 4]; //shift by 3, to be able to differentiate between 2, 3 and 4 flags nad literal run of 2, 3 or 4; 1 for sentinel
    
    uint32_t bm_comp_pos_size;
    uint32_t bm_comp_cp_size;
    uint32_t bitsize_perm;
    
    // bit vectors
    CBitMemory bm_comp_pos; //position of unique in core
    CBitMemory bm_comp_copy_orgl_id; //id of copied vector (where is original)
   
    CBitMemory bm;
    CBitMemory bv_perm;

public:
    CompressedPack()
    {
#ifdef MMAP
        fm = nullptr;
//        region = nullptr;
#endif
    }
    
    ~CompressedPack()
    {
#ifdef MMAP
        if(fm)
            delete fm;
//        if(region)
//            delete region;
#endif
        if(huf_literals)
            delete [] huf_literals;
        if(huf_zeros_runs)
            delete [] huf_zeros_runs;
        if(huf_ones_runs)
            delete [] huf_ones_runs;
        if(huf_match_lens)
            delete [] huf_match_lens;
    }
    
    bool loadPack(const std::string & arch_name);
    void getPermArray(int block_id, uint32_t * perm);
};

#endif /* compressed_pack_h */
