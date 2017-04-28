/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#ifndef decompressor_h
#define decompressor_h

#include <stdio.h>

#include "params.h"
#include "compressed_pack.h"
#include "htslib/vcf.h"
#include "samples.h"
#include "buffered_bm.h"
#include "huffman.h"
#include "my_vcf.h"
#include <deque>


class Decompressor {
    
    CompressedPack pack;
    
    string arch_name;
    
    // Samples
    Samples smpl;
    string samples_to_decompress;
    
    // Range
    string range;
    
    // Output and output settings
    htsFile *out;
    file_type out_type;
    string out_name;
    char compression_level;
    bool MB_memory;
    bool out_AC_AN;
    bool out_genotypes;
    uint32_t records_to_process;
    uint32_t max_MB_memory;
    double minAF, maxAF;
    uint32_t minAC, maxAC;
    
    // For BCF output
    htsFile * bcf = nullptr;
    hts_idx_t * bcf_idx = nullptr;
    bcf_hdr_t * hdr = nullptr;
    uint32_t * sampleIDs = nullptr;
    
    CBufferedBitMemory buff_bm;
    std::unordered_map<uint64_t, uchar_t *> done_unique;
    std::unordered_map<uint64_t, uchar_t *>::const_iterator got_it;
    
    int decompressRange(const string & range);
    void decomp_vec_rrr_range(uint64_t vec_id, uint64_t offset, uint64_t length, uint32_t & pos, uchar_t *decomp_data, uint64_t start_id, bool is_unique_id);
    
    int decompressSampleSmart(const string & range);
    
    uchar_t get_vec_byte(uint32 vec_id, uint32 byte_no, uchar_t * resUnique, bool & is_uniqe_id, uint64_t & curr_zeros, uint64_t & curr_copy);
    
    uchar_t get_vec_bytes(uint64_t vec_id, std::vector< std::pair<uint32_t, uchar_t> > & whichByte_whereInRes, uchar_t * resUnique, bool & is_uniqe_id, uint64_t & curr_zeros, uint64_t & curr_copy, uchar_t * resAll);
    
    int decompressRangeSample(const string & range);
    bool setACAN(bcf_hdr_t * hdr, bcf1_t * record, const kstring_t & str);
    
	void inline decode_perm(int no_haplotypes, int vec2_start, uint32_t *perm, uchar_t *decomp_data_perm, uchar_t *decomp_data);
	void inline decode_perm_rev(int no_haplotypes, int vec2_start, uint32_t *rev_perm, uchar_t *decomp_data_perm, uchar_t *decomp_data);
	void inline reverse_perm(uint32_t *perm, uint32_t *rev_perm, int no_haplotypes);
	
	uchar_t perm_lut[8];
    uchar lut[256][256][8];
    void initialLut();
    
    uchar_t *zeros_only_vector = nullptr;
    uchar_t *ones_only_vector = nullptr;
    int nesting = 0;
    
    // Create a deque containing ids of remembered vectors
    std::deque<uint64_t> stored_unique;
    uint64_t max_stored_unique = 0;
    
    sdsl::bit_vector zeros_only_bit_vector[2];
    sdsl::rank_support_v5<> rank_zeros_only_vector[2];
    
    sdsl::bit_vector copy_bit_vector[2];
    sdsl::rank_support_v5<> rank_copy_bit_vector[2];

    FILE * bv_out = nullptr; //for view_dev only
public:
    Decompressor()
    {
        compression_level = '1';
        out_type = VCF;
        out_name = "";
        samples_to_decompress = "";
        max_MB_memory = 0;
        MB_memory = true;
        range = "";
        records_to_process = UINT32_MAX;
        
        out_AC_AN = false;
        out_genotypes = true;
        
        minAC = 0;
        maxAC = UINT32_MAX;
        minAF = 0;
        maxAF = 1;
    }
    
    Decompressor(Params & params)
    {
        arch_name = params.arch_name;
        compression_level = params.compression_level;
        out_type = params.out_type;
        out_name = params.out_name;
        samples_to_decompress = params.samples;
        max_MB_memory = params.max_MB_memory;
        MB_memory = params.MB_memory;
        range = params.range;
        out_AC_AN = params.out_AC_AN;
        out_genotypes = params.out_genotypes;
        records_to_process = params.records_to_process;
        minAF = params.minAF;
        maxAF = params.maxAF;
        minAC = params.minAC;
        maxAC = params.maxAC;
    }
    
    ~Decompressor()
    {
        if(zeros_only_vector)
            delete [] zeros_only_vector;
        if(ones_only_vector)
            delete [] ones_only_vector;
        if(sampleIDs)
            delete [] sampleIDs;
        
        if(bcf)
            hts_close(bcf);
        if(bcf_idx)
            hts_idx_destroy(bcf_idx);
        if(hdr)
            bcf_hdr_destroy(hdr);
        
        if(bv_out)
            fclose(bv_out);
    }
    
    void decompress();
   
    
    bool loadPack();
    void setMemoryUsage();
    int initOut();
    int loadBCF();

#ifdef DEVELOPMENT_MODE
    //dev only
    int initBVOut(bool header_on);
    void getBV();
    void getVarBV(uint32_t var_to_dec);
    void getVarSampleBV(uint32_t var_to_dec, uint32_t sample_to_dec);
    void getSampleBV(uint32_t sample_to_dec);
#endif
    
};
#endif /* decompressor_h */
