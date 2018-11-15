/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include "compressed_pack.h"

using namespace std;

bool CompressedPack::loadPack(const std::string & arch_name)
{
    string fname = arch_name + ".gtc";
    
    //sdsl vectors
    sdsl::isfstream in(fname, std::ios::binary | std::ios::in);
    if (!in) 
    {
       // if (sdsl::util::verbose)
        {
            std::cerr << "Could not load file `" << fname << "`" << std::endl;
        }
        exit(1);
    }
    
    rrr_zeros_only_bit_vector[0].load(in);
    rrr_zeros_only_bit_vector[1].load(in);
    rrr_copy_bit_vector[0].load(in);
    rrr_copy_bit_vector[1].load(in);
  
    uint64_t parametersFileStartPosition = in.tellg();
    
    in.close();
    
    if (sdsl::util::verbose) {
        std::cerr << "Load file `" << fname << "`" << std::endl;
    }
    
    rrr_rank_zeros_only_bit_vector_0 = sdsl::rrr_vector<>::rank_0_type(&rrr_zeros_only_bit_vector[0]);
    rrr_rank_zeros_only_bit_vector_1 = sdsl::rrr_vector<>::rank_1_type(&rrr_zeros_only_bit_vector[1]);
    rrr_rank_copy_bit_vector[0] = sdsl::rrr_vector<>::rank_1_type(&rrr_copy_bit_vector[0]);
    rrr_rank_copy_bit_vector[1] = sdsl::rrr_vector<>::rank_1_type(&rrr_copy_bit_vector[1]);
    
//rest of archive
    uint64_t buf_pos = 0;
    uint64_t arch_size = 0;
    
#ifndef MMAP
    {
        FILE * comp = fopen(fname.c_str(), "rb");
        
        if(!comp)
        {
            cout <<"Input file ("<< fname << ")error\n";
            exit(1);
        }        
        
        fseek(comp, 0, SEEK_END);
        arch_size = ftell(comp) - parametersFileStartPosition;
        fseek(comp, parametersFileStartPosition, SEEK_SET);
        
        buf = new uchar[arch_size];
        fread(buf, sizeof(uchar)*arch_size, 1, comp);
        
        fclose(comp);
        
        cout << parametersFileStartPosition << endl;
        cout << arch_size << endl;
    }
#else //if BOOST
/*    const boost::interprocess::mode_t mode = boost::interprocess::read_only;
    fm = new boost::interprocess::file_mapping(fname.c_str(), mode);
    region = new boost::interprocess::mapped_region(*fm, mode, 0, 0);
    buf = reinterpret_cast<unsigned char*>(region->get_address()) + parametersFileStartPosition;
    arch_size = region->get_size() - parametersFileStartPosition;*/

    fm = new memory_mapped_file::read_only_mmf(fname.c_str());

    if(!fm->is_open())
    {
        cerr << "No file: " << fname << endl;
        exit(1);
    }

    buf = (uchar *) fm->data() + parametersFileStartPosition;
    arch_size = fm->file_size() - parametersFileStartPosition;
#endif
    
    memcpy(&s.ones_ranges, buf + buf_pos, sizeof(uchar));
    buf_pos = buf_pos + sizeof(uchar);

    memcpy(&s.ploidy, buf + buf_pos, sizeof(uchar));
    buf_pos = buf_pos + sizeof(uchar);
    
    memcpy(&s.vec_len, buf + buf_pos, sizeof(uint64_t));
    buf_pos = buf_pos + sizeof(uint64_t);

    s.bit_size_literal = 8;
    
    uint32_t len_huf;
    memcpy(&len_huf, buf + buf_pos, sizeof(uint32_t));
    buf_pos = buf_pos + sizeof(uint32_t);
    
    huf_match_diff_MSB.LoadTree(buf + buf_pos, len_huf);
    buf_pos = buf_pos + sizeof(uchar)*len_huf;
    
    if(huf_match_diff_MSB.max_len > 31)
    { cout << "Fast LUT impossible for huf_match_diff_MSB(huf_group_type.max_len = " << huf_match_diff_MSB.max_len << ")\n" << endl; exit(1);}
    
    memcpy(&len_huf, buf + buf_pos, sizeof(uint32_t));
    buf_pos = buf_pos + sizeof(uint32_t);
    
    huf_group_type.LoadTree(buf + buf_pos, len_huf);
    buf_pos = buf_pos + sizeof(uchar)*len_huf;
    
    if(huf_group_type.max_len > 31)
    { cout << "Fast LUT impossible for huf_group_type (huf_group_type.max_len = " << huf_group_type.max_len << ")\n" << endl; exit(1);}
    
    huf_literals = new CHuffman[s.ones_ranges];
    huf_zeros_runs = new CHuffman[s.ones_ranges];
    huf_ones_runs = new CHuffman[s.ones_ranges];
    huf_match_lens = new CHuffman[s.ones_ranges];
    for(int o_g = 0; o_g < (int) s.ones_ranges; o_g++)
    {
        memcpy(&len_huf, buf + buf_pos, sizeof(uint32_t));
        buf_pos = buf_pos + sizeof(uint32_t);
        huf_literals[o_g].LoadTree(buf + buf_pos, len_huf);
        buf_pos = buf_pos + sizeof(uchar)*len_huf;
        
        memcpy(&len_huf, buf + buf_pos, sizeof(uint32_t));
        buf_pos = buf_pos + sizeof(uint32_t);
        huf_zeros_runs[o_g].LoadTree(buf + buf_pos, len_huf);
        buf_pos = buf_pos + sizeof(uchar)*len_huf;
        
        memcpy(&len_huf, buf + buf_pos, sizeof(uint32_t));
        buf_pos = buf_pos + sizeof(uint32_t);
        huf_ones_runs[o_g].LoadTree(buf + buf_pos, len_huf);
        buf_pos = buf_pos + sizeof(uchar)*len_huf;
        
        memcpy(&len_huf, buf + buf_pos, sizeof(uint32_t));
        buf_pos = buf_pos + sizeof(uint32_t);
        huf_match_lens[o_g].LoadTree(buf + buf_pos, len_huf);
        buf_pos = buf_pos + sizeof(uchar)*len_huf;
    }
    
    // Flags for runs of literals
    memcpy(&len_huf, buf + buf_pos, sizeof(uint32_t));
    buf_pos = buf_pos + sizeof(uint32_t);
    
    huf_flags.LoadTree(buf + buf_pos, len_huf);
    buf_pos = buf_pos + sizeof(uchar)*len_huf;
    
    memcpy(&no_vec, buf + buf_pos, sizeof(uint64_t));
    buf_pos = buf_pos + sizeof(uint64_t);
    
    s.bit_size_id = (uint32_t)log2(no_vec) + 1;
    s.bit_size_match_len = (uint32_t)log2(s.vec_len) + 1;
    s.bit_size_run_len = (uint32_t)log2(s.vec_len) + 1;
    
    memcpy(&no_copy, buf + buf_pos, sizeof(uint64_t));
    buf_pos = buf_pos + sizeof(uint64_t);
    
    memcpy(&used_bits_cp, buf + buf_pos, sizeof(char));
    buf_pos = buf_pos + sizeof(char);
    
    memcpy(&bm_comp_cp_size, buf + buf_pos, sizeof(int));
    buf_pos = buf_pos + sizeof(int);
    
    bm_comp_copy_orgl_id.Open(buf + buf_pos, bm_comp_cp_size);
    buf_pos = buf_pos + sizeof(uchar_t) * bm_comp_cp_size;
    
    memcpy(&no_non_copy, buf + buf_pos, sizeof(uint64_t));
    buf_pos = buf_pos + sizeof(uint64_t);
    
    for(int o = 0; o < (int) s.ones_ranges; o++)
    {
        memcpy(&max_used_bits_litRunSize[o], buf + buf_pos, sizeof(uchar_t));
        buf_pos = buf_pos + sizeof(uchar_t);
        
        memcpy(&used_bits_litRunSize[o], buf + buf_pos, sizeof(uchar_t));
        buf_pos = buf_pos + sizeof(uchar_t);
        
        memcpy(minLitRunSize[o], buf + buf_pos, sizeof(uint32_t) * (MAX_LITERAL_RUN + 4));
        buf_pos = buf_pos + sizeof(uint32_t) * (MAX_LITERAL_RUN + 4);
    }
    
    memcpy(&used_bits_noncp, buf + buf_pos, sizeof(char));
    buf_pos = buf_pos + sizeof(char);
    
    memcpy(&s.bit_size_id_match_pos_diff, buf + buf_pos, sizeof(char));
    buf_pos = buf_pos + sizeof(char);
    
    memcpy(&bm_comp_pos_size, buf + buf_pos, sizeof(bm_comp_pos.mem_buffer_pos));
    buf_pos = buf_pos + sizeof(bm_comp_pos.mem_buffer_pos);
    
    bm_comp_pos.Open( buf + buf_pos, bm_comp_pos_size);
    buf_pos = buf_pos + sizeof(uchar_t)*bm_comp_pos_size;
    
    memcpy(&no_blocks, buf + buf_pos, sizeof(uint32_t));
    buf_pos = buf_pos + sizeof(uint32_t);

    memcpy(&s.max_no_vec_in_block, buf + buf_pos, sizeof(uint32_t));
    buf_pos = buf_pos + sizeof(uint32_t);

    memcpy(&s.n_samples, buf + buf_pos, sizeof(uint32_t));
    buf_pos = buf_pos + sizeof(uint32_t);

    uint64_t bv_perm_size;
    memcpy(&bv_perm_size, buf + buf_pos, sizeof(bv_perm.mem_buffer_pos));
    buf_pos = buf_pos + sizeof(bv_perm.mem_buffer_pos);
    bv_perm.Open(buf + buf_pos, bv_perm_size);
    buf_pos += bv_perm_size;
    
    uint64_t core_size = arch_size - buf_pos;
    bm.Open(buf + buf_pos, core_size);

    return true;
}

void CompressedPack::getPermArray(int block_id, uint32_t * perm)
{
    uint32_t no_haplotypes = s.n_samples * s.ploidy;
    uint32_t bits_used_single = s.bits_used(no_haplotypes);
    uint32_t single_perm_bv_size = bits_used_single * no_haplotypes; //in bits
    single_perm_bv_size = single_perm_bv_size/8 + (single_perm_bv_size%8?1:0); //in bytes
    bv_perm.SetPos(block_id * single_perm_bv_size);
    
    for(uint32_t i = 0; i < no_haplotypes; ++i)
    {
        if(!bv_perm.GetBits(perm[i], bits_used_single))
        {
            cout << "error in getPermArray" << endl;
            exit(1);
        }
    }    
}
