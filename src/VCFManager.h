/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#ifndef VCFManager_h
#define VCFManager_h

#include <iostream>
#include <fstream>
#include "htslib/vcf.h"
#include "defs.h"
#include "params.h"
#include "bit_memory.h"
#include "queues.h"

class VCFManager {
    
    htsFile *in_vcf = nullptr;
    htsFile *out_vcf = nullptr;
    
    bool in_open, out_open;
    file_type in_type, out_type;
    std::string in_vcf_name, out_vcf_name, arch_name;
    
    bcf_hdr_t *in_hdr = nullptr, *out_hdr = nullptr;

    uint32_t no_samples;
    uint32_t ploidy;
    uint64_t vec_len;
    
    uint64_t no_vec;
    
    CBitMemory bv;
    int64 block_max_size;
    uint32_t no_vec_in_block, vec_read_in_block, block_id;
    CBlockQueue * queue = nullptr;
    
    bool in_hdr_read;
    
    void ReadInHdr();
    bool OpenInVCF();
    bool OpenOutVCF();
    void setBitVector();
    void addGTtoBitVector(bcf_hdr_t * hdr, bcf1_t * rec);
    
public:
    VCFManager() {
        in_open = false;
        out_open = false;
        in_hdr_read = false;
        no_samples = 0;
        ploidy = 0;
        no_vec = 0;
    }
    VCFManager(const Params & params) {
        in_open = false;
        out_open = false;
        in_hdr_read = false;
        no_samples = 0;
        ploidy = params.ploidy;
        
        in_vcf_name = params.in_file_name;
        in_type = params.in_type;
        
        arch_name = params.arch_name;
        no_vec = 0;
        
        if(params.task == tcompress || params.task == tcompress_dev_pre)
        {
            out_vcf_name = arch_name;
            out_vcf_name += ".bcf";
            out_type = BCF;
        }
        else if(params.task == tquery || params.task == tquery_dev)
        {
            out_type = params.out_type;
            if(out_type == VCF)
                out_vcf_name = params.out_name + ".vcf";
            else
                out_vcf_name = params.out_name + ".bcf";
        }
        
        if(params.task == tcompress || params.task == tcompress_dev)
        {
            no_vec_in_block = params.var_in_block*2;
        }
    }
    
    void CloseFiles()
    {
        if(in_open)
        {
            bcf_hdr_destroy(in_hdr);
            in_hdr = nullptr;
            hts_close(in_vcf);
            in_vcf = nullptr;
            in_open = false;
        }
        if(out_open)
        {
            
            bcf_hdr_destroy(out_hdr);
            out_hdr = nullptr;
            hts_close(out_vcf);
            out_vcf = nullptr;
            out_open = false;
        }
    }
    
    ~VCFManager() {
        CloseFiles();
    }
  
    void setQueue(CBlockQueue * _queue)
    {
        queue = _queue;
    };
    
    uint64_t getNoVec()
    {
        return no_vec;
    }
    
    uint32_t getInputNoSamples();
    bool CreateIndividualListFile();
    bool ProcessInVCF();
    
};

#endif /* VCFreader_hpp */
