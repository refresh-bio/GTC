/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#ifndef end_compressor_h
#define end_compressor_h

#include <iostream>
#include <vector>
#include "params.h"
#include "defs.h"
#include "compression_settings.h"
#include <sdsl/bit_vectors.hpp>
#include "huffman.h"
#include "bit_memory.h"

//#define MMAP

//#ifdef MMAP
//#include <boost/interprocess/file_mapping.hpp>
//#include <boost/interprocess/mapped_region.hpp>
//#include <cpp-mmf/memory_mapped_file.hpp>
//#endif

class EndCompressor{
    
    CompSettings * s = nullptr;
    uint64_t no_vec, curr_vec_id = 0;
    sdsl::bit_vector zeros_only_bit_vector[2];
    sdsl::bit_vector copy_bit_vector[2];
    sdsl::bit_vector unique;
    sdsl::rank_support_v5<>  rank_unique;
    sdsl::rank_support_v5<> rank_zeros_only_vector[2];
    sdsl::rank_support_v5<> rank_copy_bit_vector[2];
    
    uint32_t curr_pos = 0;
    uint64_t copy_no = 0, unique_no = 0;
    vector<uint64_t> unique_in_block;
    
    uint32_t * comp_pos_non_copy;
    uint32_t * comp_pos_copy;
    
    CBitMemory bm;  //initial
    CBitMemory bm_huff;
    CBitMemory bm_comp_pos;
    CBitMemory bm_comp_copy_orgl_id;
    
    CHuffman * huf_literals = NULL;
    CHuffman * huf_zeros_runs = NULL;
    CHuffman * huf_ones_runs = NULL;
    
    CHuffman * huf_match_lens = NULL;
    CHuffman huf_flags;
    CHuffman huf_group_type;
    CHuffman huf_match_diff_MSB;

    
    uint64_t * hist_match_diff_MSB = nullptr; //
    uint64_t ** hist_literals  = NULL;//
    uint64_t * hist_group_type  = NULL;//
    uint64_t ** hist_zero_runs = NULL;//
    uint64_t ** hist_ones_runs = NULL;//
    
    uint64_t ** hist_match_lens = NULL;//
    uint64_t hist_flags[MAX_LITERAL_RUN+4] = { 0 }; //shift by 3, to be able to difference
    uint32_t hist_litRunSize_usedBits[MAX_NUMBER_OF_GROUP][32]= { {0} };
   // uint32_t hist_litRunSize[MAX_NUMBER_OF_GROUP][2048]= { {0} };
   
    uint32_t literalRunCount = 0;
    uint16_t * litRunSize = NULL;
    uint32_t minLitRunSize[MAX_NUMBER_OF_GROUP][MAX_LITERAL_RUN+4]; //shift by 3, to be able to differentiate between 2, 3 and 4 flags nad literal run of 2, 3 or 4
    uint32_t  max_used_bits_litRunSize[MAX_NUMBER_OF_GROUP], used_bits_litRunSize[MAX_NUMBER_OF_GROUP] ;    uint32_t litRunSizeMax[MAX_NUMBER_OF_GROUP] = { 0 };

    uint32_t used_bits_cp;
    uint32_t  used_bits_noncp;
     uint32_t max_pos_diff = 0;
    uint32_t  no_blocks = 0;
    uint32_t block = 0;
    
    std::vector< std::vector<int> > perms;
    
    void encode_vec(uint64_t vec_id);
    
    char bits_used(unsigned int n) ;
    
    void calcModel();
    
    void getLitRunSizes(uint64_t vec_id);
    void decreaseLitRunSizesbyMin(uint64_t unique_no);
    
public:
    int storeArchive(const char * arch_name);
    
    ~EndCompressor()
    {
        if(hist_match_diff_MSB)
            delete [] hist_match_diff_MSB;
        
        if(hist_group_type)
            delete [] hist_group_type;
        if(hist_zero_runs)
        {
            for(int o = 0; o < (int) s->ones_ranges; o++)
            {
                if(hist_zero_runs[o])
                    delete [] hist_zero_runs[o];
                
            }
            delete [] hist_zero_runs;
        }
        
        if(hist_ones_runs)
        {
            for(int o = 0; o < (int) s->ones_ranges; o++)
            {
                if(hist_ones_runs[o])
                    delete [] hist_ones_runs[o];
                
            }
            delete [] hist_ones_runs;
        }
        
        
        if(hist_literals)
        {
            for(int o = 0; o < (int) s->ones_ranges; o++)
            {
                if(hist_literals[o])
                    delete [] hist_literals[o];
            }
            delete [] hist_literals;
        }
        
        if(hist_match_lens)
        {
            for(int o = 0; o < (int) s->ones_ranges; o++)
            {
                if(hist_match_lens[o])
                    delete [] hist_match_lens[o];
            }
            delete [] hist_match_lens;
        }
        
//#ifndef MMAP
//        if(buf)
//            delete [] buf;
//#endif
        
        if(litRunSize)
            delete [] litRunSize;
        
        if(huf_literals)
            delete [] huf_literals;
        if(huf_match_lens)
            delete [] huf_match_lens;
        if(huf_zeros_runs)
            delete [] huf_zeros_runs;
        if(huf_ones_runs)
            delete [] huf_ones_runs;
        
        
        huf_flags.Restart();
        huf_group_type.Restart();
        huf_match_diff_MSB.Restart();
    }
    
    EndCompressor(CompSettings * _settings, uint64_t _no_vec)
    {
        s = _settings;
        no_vec = _no_vec;
        zeros_only_bit_vector[0] = sdsl::bit_vector(no_vec/2+no_vec%2, 0);
        zeros_only_bit_vector[1] = sdsl::bit_vector(no_vec/2+no_vec%2, 0);
        copy_bit_vector[0] = sdsl::bit_vector(no_vec/2+no_vec%2, 0);
        copy_bit_vector[1] = sdsl::bit_vector(no_vec/2+no_vec%2, 0);
        
        curr_vec_id = 0;
        curr_pos = 0;
        copy_no = 0;
        unique_no = 0;
        no_blocks = 0;

        
        comp_pos_non_copy = new uint32_t[no_vec]();
        comp_pos_copy = new uint32_t[no_vec]();
        
        bm.Create((no_vec*s->n_samples*s->ploidy/8)/10);
        
        unique = sdsl::bit_vector(no_vec, 0);
        
        perms.clear();
    }
    void AddBlock(int &id_block, unsigned char *compressed_block, size_t n_recs, size_t compressed_size, std::vector<int> &perm, std::vector<bool> &zeros, std::vector<bool> &copies, uint32_t * origin_of_copy);
    void Encode();
    
};

#endif /* end_compressor_h */
