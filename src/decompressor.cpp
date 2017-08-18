/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include "decompressor.h"

void Decompressor::decompress()
{
    if(samples_to_decompress == "")
        decompressRange(range);
    else
        decompressSampleSmart(range);
}

int Decompressor::decompressRange(const string & range)
{
    initialLut();
    
    uchar_t * decomp_data = nullptr;
    uchar_t * decomp_data_perm = nullptr;
    uint32_t * perm = nullptr;
    uint32_t * rev_perm = nullptr;
    uint32_t pos;
    
    uint32_t no_haplotypes = smpl.no_samples * pack.s.ploidy;
    perm = new uint32_t[no_haplotypes];
    rev_perm = new uint32_t[no_haplotypes];
	
    bcf1_t * record = bcf_init();
    uint32_t g, vec1_start, vec2_start;
    
    uint32_t  end;
    uint32_t i = 0;
    
    buff_bm.setBitMemory(&pack.bm);
    
    const char *key = "GT";
    khint_t k;
    vdict_t *d;
    int fmt_id;
    uint32_t block_id, prev_block_id = 0xFFFF;
    done_unique.clear();
    kstring_t str = {0,0,0};
    uint32_t written_records = 0;    
   
    long long * tmp_vec_ll = nullptr;
    long long * lut_ll = nullptr;
    if(range == "")
    {
        // VCF/BCF processing
        d = (vdict_t*)hdr->dict[BCF_DT_ID];
        k = kh_get(vdict, d, key);
        fmt_id = (k == kh_end(d)? -1 : kh_val(d, k).id);
        
        if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,fmt_id) )
        {
            if(perm)
                delete [] perm;
            if(rev_perm)
                delete [] rev_perm;
            bcf_destroy(record);
            
            if ( !no_haplotypes )
                return 0;
            return -1;  // the key not present in the header
        }
        
        bcf_enc_int1(&str, fmt_id);
        bcf_enc_size(&str, pack.s.ploidy, BCF_BT_INT8);
    
        char *tmp;
        str.m = str.l + no_haplotypes + 1;
        kroundup32(str.m);
        if ((tmp = (char*)realloc(str.s, str.m)))
            str.s = tmp;
        else
            exit(8);
        
        str.l = 3;
        if (no_haplotypes == 0) bcf_enc_size(&str, 0, BCF_BT_NULL);
        
        if(no_haplotypes & 7)//%8)
        {
            end = pack.s.vec_len - 1;
            g = no_haplotypes & 7;//%8;
        }
        else
        {
            end = pack.s.vec_len;
            g = 0;
        }
        
        decomp_data = new uchar_t[pack.s.vec_len*2];
        decomp_data_perm = new uchar_t[pack.s.vec_len*2];
        
        while (bcf_read1(bcf, hdr, record)>= 0 && written_records < records_to_process)
        {
            pos = 0;
            str.l = 3;
            
            bcf_unpack(record, BCF_UN_ALL);
            record->n_sample = bcf_hdr_nsamples(hdr);
            
            // Permutations
            block_id = i/pack.s.max_no_vec_in_block; // Pair of vectors always in the same block (so it is enough to find first vector block
            if(block_id != prev_block_id)
            {
                // Get approproate permutations
                pack.getPermArray(block_id, perm);
                reverse_perm(perm, rev_perm, no_haplotypes);
            }
            
            fill_n(decomp_data, pack.s.vec_len*2, 0);
            
            vec1_start = 0;
            vec2_start = (uint32_t)pack.s.vec_len;
            decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
            decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
            
            decode_perm_rev(no_haplotypes, vec2_start, rev_perm, decomp_data_perm, decomp_data);
        
            prev_block_id = block_id;
           
            char *pt = str.s + str.l;
            tmp_vec_ll = (long long *)pt;
            for (vec1_start = 0; vec1_start < end; ++vec1_start)
            {
                //memcpy(pt + (vec1_start << 3), lut[decomp_data[vec1_start]][decomp_data[vec2_start++]], 8);
                lut_ll = (long long *)(lut[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
               *(tmp_vec_ll+ vec1_start) = *lut_ll;
            }
            
            str.l = str.l + (end << 3);
            if(g)
            {
                memcpy(str.s + str.l, lut[decomp_data[vec1_start]][decomp_data[vec2_start]], g);
                str.l = str.l + g;
            }
            
            str.s[str.l] = 0;
            
            bcf_update_info_int32(hdr, record, "_row", NULL, 0);
            
            // AC/AN count
            if(out_AC_AN)
            {
                if(setACAN(hdr, record, str))
                {
                    if(out_genotypes)
                        bcf_update_genotypes_fast(str, record);
                    
                    bcf_write1(out, hdr, record);
                    written_records++;
                }
            }
            else
            {
                if(out_genotypes)
                    bcf_update_genotypes_fast(str, record);
            
                bcf_write1(out, hdr, record);
                written_records++;
            }            
        }        
        
        for (auto & it:done_unique)
        {
            delete [] it.second;
        }

        done_unique.clear();
    }
    else
    {
        // VCF/BCF processing
        d = (vdict_t*)hdr->dict[BCF_DT_ID];
        k = kh_get(vdict, d, key);
        fmt_id = (k == kh_end(d)? -1 : kh_val(d, k).id);
        
        if (!bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,fmt_id))
        {
            if(perm)
                delete [] perm;
            if(rev_perm)
                delete [] rev_perm;
            bcf_destroy(record);
            
            if ( !no_haplotypes ) return 0;
            return -1;  // the key not present in the header
        }
              
        bcf_enc_int1(&str, fmt_id);
        bcf_enc_size(&str, pack.s.ploidy, BCF_BT_INT8);
 
        char *tmp;
        str.m = str.l + no_haplotypes + 1;
        kroundup32(str.m);
        if ((tmp = (char*)realloc(str.s, str.m)))
            str.s = tmp;
        else
            exit(8);
    
        str.l = 3;
        if (no_haplotypes == 0) bcf_enc_size(&str, 0, BCF_BT_NULL);
        
        if(no_haplotypes & 7)//%8)
        {
            end = pack.s.vec_len - 1;
            g = no_haplotypes & 7;//%8;
        }
        else
        {
            end = pack.s.vec_len;
            g = 0;
        }
        
        decomp_data = new uchar_t[pack.s.vec_len*2];
        decomp_data_perm = new uchar_t[pack.s.vec_len*2];
       
        bcf_info_t * a;
        hts_itr_t * itr = bcf_itr_querys(bcf_idx, hdr, range.c_str());

        while(bcf_itr_next(bcf, itr, record) != -1 && written_records < records_to_process)
        {
            if(!i)
            {
                a = bcf_get_info(hdr, record, "_row"); //bcf_get_info_id
                i = (a->v1.i*2);
            }
            
            pos = 0;
            str.l = 3;
            
            bcf_unpack(record, BCF_UN_ALL);
            record->n_sample = bcf_hdr_nsamples(hdr);

            // Permutations
            block_id = i/pack.s.max_no_vec_in_block; // Pair of vectors always in the same block (so it is enough to find first vector block
            if(block_id != prev_block_id)
            {
                // Get approproate permutations
                pack.getPermArray(block_id, perm);
                reverse_perm(perm, rev_perm, no_haplotypes);
            }
            fill_n(decomp_data, pack.s.vec_len*2, 0);

            vec1_start = 0;
            vec2_start = pack.s.vec_len;
            decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
            decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
            
            // Permutations
            decode_perm_rev(no_haplotypes, vec2_start, rev_perm, decomp_data_perm, decomp_data);
            
            prev_block_id = block_id;
            char *pt = str.s + str.l;
            tmp_vec_ll = (long long *)pt;
            for (vec1_start = 0; vec1_start < end; ++vec1_start)
            {
                //memcpy(pt + (vec1_start << 3), lut[decomp_data[vec1_start]][decomp_data[vec2_start++]], 8);
                lut_ll = (long long *)(lut[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                *(tmp_vec_ll+ vec1_start) = *lut_ll;
            }
            
            str.l = str.l + (end << 3);
            if(g)
            {
                memcpy(str.s + str.l, lut[decomp_data[vec1_start]][decomp_data[vec2_start]], g);
                str.l = str.l + g;
            }
            
            str.s[str.l] = 0;

            bcf_update_info_int32(hdr, record, "_row", NULL, 0);
       
            // AC/AN count
            if(out_AC_AN)
            {
                if(setACAN(hdr, record, str))
                {
                    if(out_genotypes)
                        bcf_update_genotypes_fast(str, record);
                    bcf_write1(out, hdr, record);
                    written_records++;
                }
            }
            else
            {
                if(out_genotypes)
                    bcf_update_genotypes_fast(str, record);
                bcf_write1(out, hdr, record);
                written_records++;
            }
            
        }
        bcf_itr_destroy(itr);
        
        for (auto & it:done_unique)
        {
            delete [] it.second;
        }
        done_unique.clear();
    }
    
    hts_close(out);
    
    if(decomp_data)
        delete [] decomp_data;
    if(decomp_data_perm)
        delete [] decomp_data_perm;
    if(perm)
        delete [] perm;
    if(rev_perm)
        delete [] rev_perm;
    if(str.s)
    {
        free(str.s);
        str.s = nullptr;
    }
    bcf_destroy(record);
    
    return 0;
}

int Decompressor::decompressSampleSmart(const string & range)
{
    if(smpl.no_samples > NO_SAMPLE_THRESHOLD)// || range != "")
    {
        return decompressRangeSample(range);
    }
    initialLut();
    
    std::vector<std::pair<uint32_t,uint32_t>> whichByte_whereInRes(smpl.no_samples*pack.s.ploidy + 1); // +1 for guard
    //std::vector<std::pair<uint32_t,uchar_t>> whichByte_whereInRes(smpl.no_samples*pack.s.ploidy);
    
    buff_bm.setBitMemory(&pack.bm);
    bool is_unique = false;
    uint32_t unique_pos = 0, unique_pos_first_in_block = 0;
    uint64_t  curr_zeros = 0, curr_copy = 0;
    
    uint64_t bv_size = pack.rrr_zeros_only_bit_vector[0].size();
    
    zeros_only_bit_vector[0] = sdsl::bit_vector(bv_size);
    zeros_only_bit_vector[1] = sdsl::bit_vector(bv_size);
    copy_bit_vector[0] = sdsl::bit_vector(bv_size);
    copy_bit_vector[1] = sdsl::bit_vector(bv_size);
    
    uint64_t v_pos;
    for(v_pos = 0; v_pos + 64 < bv_size; v_pos += 64)
    {
        zeros_only_bit_vector[0].set_int(v_pos, pack.rrr_zeros_only_bit_vector[0].get_int(v_pos, 64), 64);
        zeros_only_bit_vector[1].set_int(v_pos, pack.rrr_zeros_only_bit_vector[1].get_int(v_pos, 64), 64);
        copy_bit_vector[0].set_int(v_pos, pack.rrr_copy_bit_vector[0].get_int(v_pos, 64), 64);
        copy_bit_vector[1].set_int(v_pos, pack.rrr_copy_bit_vector[1].get_int(v_pos, 64), 64);
    }
    
    uint64_t tail_len = bv_size - v_pos;
    zeros_only_bit_vector[0].set_int(v_pos, pack.rrr_zeros_only_bit_vector[0].get_int(v_pos, tail_len), tail_len);
    zeros_only_bit_vector[1].set_int(v_pos, pack.rrr_zeros_only_bit_vector[1].get_int(v_pos, tail_len), tail_len);
    copy_bit_vector[0].set_int(v_pos, pack.rrr_copy_bit_vector[0].get_int(v_pos, tail_len), tail_len);
    copy_bit_vector[1].set_int(v_pos, pack.rrr_copy_bit_vector[1].get_int(v_pos, tail_len), tail_len);
    
    uint32_t * perm = nullptr;
    
    perm = new uint32_t[pack.s.n_samples * pack.s.ploidy]; //pack.s.n_samples * pack.s.ploidy
    
    uint32_t block_id, prev_block_id = 0xFFFF;

    uint32_t where; //uchar where;
    uint32_t ind_id_orig;
    
    unique_pos = 0;
    curr_zeros = 0, curr_copy = 0;
  
    int32_t *tmpia = (int*)malloc(smpl.no_samples * pack.s.ploidy*sizeof(int));
    bcf1_t * record = bcf_init();
    uint64_t var_idx = 0, start;
    uint32_t written_records = 0;
    
    // Set out BCF
    const char *key = "GT";
    khint_t k;
    vdict_t *d;
    int fmt_id;
  
    d = (vdict_t*)hdr->dict[BCF_DT_ID];
    k = kh_get(vdict, d, key);
    fmt_id = (k == kh_end(d)? -1 : kh_val(d, k).id);
    
    uint32_t no_haplotypes = smpl.no_samples* pack.s.ploidy;
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,fmt_id) )
    {
        hts_close(out);
        free(tmpia);
        if ( !no_haplotypes ) return 0;
        return -1;  // the key not present in the header
    }
    
    kstring_t str = {0,0,0};
    bcf_enc_int1(&str, fmt_id);
    
    bcf_enc_size(&str, pack.s.ploidy, BCF_BT_INT8);
    
    char *tmp;
    str.m = str.l + no_haplotypes + 1;
    kroundup32(str.m);
    if ((tmp = (char*)realloc(str.s, str.m)))
        str.s = tmp;
    else
        exit(8);

    str.l = 3;
//    uint32_t pos;
    if (no_haplotypes == 0) bcf_enc_size(&str, 0, BCF_BT_NULL);
    
    uchar_t *resUnique = nullptr;
   // resUnique = new uchar_t[pack.no_non_copy*smpl.no_samples*pack.s.ploidy];
    resUnique = new uchar_t[pack.s.max_no_vec_in_block*smpl.no_samples*pack.s.ploidy];
    uchar_t * resAll = nullptr;
    // resAll = new uchar_t[pack.no_vec*smpl.no_samples*pack.s.ploidy];
  //  resAll = new uchar_t[pack.s.max_no_vec_in_block*smpl.no_samples*pack.s.ploidy];
    resAll = new uchar_t[2*smpl.no_samples*pack.s.ploidy];

    
    uint32_t first_vec_in_block = 0;
    
    if(range != "")
    {
        bcf_info_t * a = nullptr;
        hts_itr_t * itr = bcf_itr_querys(bcf_idx, hdr, range.c_str());
        
        while(bcf_itr_next(bcf, itr, record) != -1 && written_records < records_to_process)
        {
            if(!var_idx)
            {
                a = bcf_get_info(hdr, record, "_row"); //bcf_get_info_id
                var_idx = a->v1.i;
            }
            
            // Permutations
            block_id = (var_idx*2)/pack.s.max_no_vec_in_block;
            
            if(block_id != prev_block_id) // Get perm and find out which bytes need decoding
            {
                first_vec_in_block = block_id*pack.s.max_no_vec_in_block;
                
                unique_pos = 0;
                if(prev_block_id == 0xFFFF) //to set curr_zeros, curr_copy (first processed block)
                {
                    uint8_t parity = first_vec_in_block&1;
                    uint32_t id = first_vec_in_block>>1;
                    curr_zeros = pack.rrr_rank_zeros_only_bit_vector_0(id+((parity))) + pack.rrr_rank_zeros_only_bit_vector_1(id);
                    curr_copy = pack.rrr_rank_copy_bit_vector[0](id+((parity))) + pack.rrr_rank_copy_bit_vector[1](id);
                    //unique_pos_first_in_block = first_vec_in_block - curr_zeros - curr_copy;
                    
                   // first_vec_in_block = 0;
                }
                unique_pos_first_in_block = first_vec_in_block - curr_zeros - curr_copy;
                
                // Get approproate permutations
                pack.getPermArray(block_id, perm);
                
                for(uint32_t s = 0; s < smpl.no_samples; s++)
                {
                    for (uint32_t p = 0; p < pack.s.ploidy; p++ )
                    {
                        ind_id_orig = sampleIDs[s]*pack.s.ploidy+p;
                        
                        where = s*pack.s.ploidy  + p;
                        whichByte_whereInRes[where] =  std::make_pair(perm[ind_id_orig]>>3, where); //now original_id_only
                    }
                }
                whichByte_whereInRes[smpl.no_samples*pack.s.ploidy] =  std::make_pair(0xFFFFFFFF, smpl.no_samples*pack.s.ploidy); //guard
                
                // Sort by byte_no
                std::sort(whichByte_whereInRes.begin(), whichByte_whereInRes.end());
                
            /*    //count bytes to decode
                uint32_t bytesToDecode = 0;
                if(no_haplotypes > 0)
                    bytesToDecode = 1;
                for(int b = 1; b < no_haplotypes; b++)
                    if(whichByte_whereInRes[b].first != whichByte_whereInRes[b-1].first)
                        bytesToDecode++;
                
              */
                prev_block_id = block_id;
                
                // Get vectors from all block
                //copies only within the same block, so only part of resUnique (within the same block and using the same perm) will be used - no need to fix previous unique bytes (these are in different blocks, with diff perm)
                uint64_t last_vec_in_block = (block_id+1)* pack.s.max_no_vec_in_block ;
                last_vec_in_block = last_vec_in_block < pack.no_vec ? last_vec_in_block : pack.no_vec;
                
                for (uint64_t i = first_vec_in_block; i < var_idx*2; i++ )
                {
                    get_vec_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy, resAll, unique_pos_first_in_block, first_vec_in_block, false);
                
                    if(is_unique)
                    {
                      //  for(uint32_t idx = 0; idx < no_haplotypes; idx++)
                        {
                           // resUnique[pack.no_non_copy*whichByte_whereInRes[idx].second  + unique_pos] = resAll[whichByte_whereInRes[idx].second * pack.no_vec + i];
                            
                       //     resUnique[unique_pos*pack.s.vec_len + whichByte_whereInRes[idx].second] = resAll[i * pack.s.vec_len + whichByte_whereInRes[idx].second];
                       //      memcpy(resUnique + unique_pos*no_haplotypes, resAll + (i-first_vec_in_block) * no_haplotypes, no_haplotypes);
                                  memcpy(resUnique + unique_pos*no_haplotypes, resAll + (i&1) * no_haplotypes, no_haplotypes);
                            
                        }
                        unique_pos++;
                    }
                }
            }
         //   else
            {
                
                //previous vectors in block are decompressed already
                for (uint64_t i = var_idx*2; i <= var_idx*2+1; i++ )
                    
                {    get_vec_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy, resAll, unique_pos_first_in_block, first_vec_in_block);
                
                    if(is_unique)
                    {
                        //  for(uint32_t idx = 0; idx < no_haplotypes; idx++)
                        {
                            // resUnique[pack.no_non_copy*whichByte_whereInRes[idx].second  + unique_pos] = resAll[whichByte_whereInRes[idx].second * pack.no_vec + i];
                            
                            //     resUnique[unique_pos*pack.s.vec_len + whichByte_whereInRes[idx].second] = resAll[i * pack.s.vec_len + whichByte_whereInRes[idx].second];
                           // memcpy(resUnique + unique_pos*no_haplotypes, resAll + (i-first_vec_in_block) * no_haplotypes, no_haplotypes);
                             memcpy(resUnique + unique_pos*no_haplotypes, resAll + (i&1) * no_haplotypes, no_haplotypes);
                        }
                        unique_pos++;
                    }
                }
                
            }
            
//            pos = 0;
            str.l = 3;
            bcf_unpack(record, BCF_UN_ALL);
            record->n_sample = bcf_hdr_nsamples(hdr);
            char *pt = str.s + str.l;
    
            for (int g = 0; g < (int) smpl.no_samples; g++)
            {
                for (uint32_t p = 0; p < pack.s.ploidy; p++ )
                {
                    start =  (g*pack.s.ploidy+p);
                    pt[start] = lut[resAll[start]][resAll[start+no_haplotypes]][(perm[sampleIDs[g]*pack.s.ploidy + p])%8];
                    //memcpy(pt + (g*pack.s.ploidy+p), lut[resAll[start]][resAll[start+no_haplotypes]] +(perm[sampleIDs[g]*pack.s.ploidy + p])%8, 1);
                }
            }
          
            str.l = str.l + (smpl.no_samples*pack.s.ploidy);
            str.s[str.l] = 0;
            
            bcf_update_info_int32(hdr, record, "_row", NULL, 0);
            
            // AC/AN count
            if(out_AC_AN)
            {
                if(setACAN(hdr, record, str))
                {
                    if(out_genotypes)
                        bcf_update_genotypes_fast(str, record);
                    
                    bcf_write1(out, hdr, record);
                    written_records++;
                    
                }
            }
            else
            {
                if(out_genotypes)
                    bcf_update_genotypes_fast(str, record);
                
                bcf_write1(out, hdr, record);
                written_records++;
            }
            
            var_idx++;
            
            prev_block_id = block_id;
        }
        bcf_itr_destroy(itr);
    }
    else
    {
        
        while (bcf_read1(bcf, hdr, record) >= 0 &&  written_records < records_to_process)
        {
            // Permutations
            block_id = (var_idx*2) / pack.s.max_no_vec_in_block;
            if(block_id != prev_block_id) //get perm and find out which bytes need decoding
            {
                first_vec_in_block = block_id*pack.s.max_no_vec_in_block;
                unique_pos = 0;
                if(prev_block_id == 0xFFFF) //to set curr_zeros, curr_copy (first processed block)
                {
                    uint8_t parity = first_vec_in_block&1;
                    uint32_t id = first_vec_in_block>>1;
                    curr_zeros = pack.rrr_rank_zeros_only_bit_vector_0(id+((parity))) + pack.rrr_rank_zeros_only_bit_vector_1(id);
                    curr_copy = pack.rrr_rank_copy_bit_vector[0](id+((parity))) + pack.rrr_rank_copy_bit_vector[1](id);
                    //unique_pos_first_in_block = first_vec_in_block - curr_zeros - curr_copy;
                    
                    // first_vec_in_block = 0;
                }
                unique_pos_first_in_block = first_vec_in_block - curr_zeros - curr_copy;
                
                
                
                //get approproate permutations
                pack.getPermArray(block_id, perm);
                
                for(uint32_t s = 0; s < smpl.no_samples; s++)
                {
                    for (uint32_t p = 0; p < pack.s.ploidy; p++ )
                    {
                        ind_id_orig = sampleIDs[s]*pack.s.ploidy+p;
                        where = s*pack.s.ploidy  + p;
                        whichByte_whereInRes[where] =  std::make_pair(perm[ind_id_orig]>>3, where); // now original_id_only
                    }
                }
                 whichByte_whereInRes[smpl.no_samples*pack.s.ploidy] =  std::make_pair(0xFFFFFFFF, smpl.no_samples*pack.s.ploidy); //guard
                
                
                // Sort by byte_id (.first)
                std::sort(whichByte_whereInRes.begin(), whichByte_whereInRes.end());
                
            /*    //count bytes to decode
                uint32_t bytesToDecode = 0;
                if(no_haplotypes > 0)
                    bytesToDecode = 1;
                for(int b = 1; b < no_haplotypes; b++)
                    if(whichByte_whereInRes[b].first != whichByte_whereInRes[b-1].first)
                        bytesToDecode++;
              */
                prev_block_id = block_id;
                
                // Get vectors from all block
                uint64_t last_vec_in_block = var_idx*2 + pack.s.max_no_vec_in_block;
                last_vec_in_block = last_vec_in_block < pack.no_vec ? last_vec_in_block : pack.no_vec;
                
                for (uint64_t i = var_idx*2; i <= var_idx*2+1; i++ ){
                    get_vec_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy, resAll, unique_pos_first_in_block, first_vec_in_block);
                    
                    
                    if(is_unique)
                    {
                      //  for(uint32_t idx = 0; idx < no_haplotypes; idx++)
                        {
                        //    resUnique[pack.no_non_copy*whichByte_whereInRes[idx].second  + unique_pos] = resAll[whichByte_whereInRes[idx].second * pack.no_vec + i];
                       //    memcpy(resUnique + unique_pos*no_haplotypes, resAll + (i-first_vec_in_block) * no_haplotypes, no_haplotypes);
                            memcpy(resUnique + unique_pos*no_haplotypes, resAll + (i&1) * no_haplotypes, no_haplotypes);

                            
                        }
                        unique_pos++;
                    }
                }
            }
            else
            {
                for (uint64_t i = var_idx*2; i <= var_idx*2+1; i++ ){
                    get_vec_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy, resAll, unique_pos_first_in_block, first_vec_in_block);
                    
                    
                    if(is_unique)
                    {
                        //  for(uint32_t idx = 0; idx < no_haplotypes; idx++)
                        {
                            //    resUnique[pack.no_non_copy*whichByte_whereInRes[idx].second  + unique_pos] = resAll[whichByte_whereInRes[idx].second * pack.no_vec + i];
                            //    memcpy(resUnique + unique_pos*no_haplotypes, resAll + (i-first_vec_in_block) * no_haplotypes, no_haplotypes);
                            memcpy(resUnique + unique_pos*no_haplotypes, resAll + (i&1) * no_haplotypes, no_haplotypes);
                            
                            
                        }
                        unique_pos++;
                    }
                }
            }
            
//            pos = 0;
            str.l = 3;
            bcf_unpack(record, BCF_UN_ALL);
            record->n_sample = bcf_hdr_nsamples(hdr);
            
            char *pt = str.s + str.l;
            
            for (int g = 0; g < (int) smpl.no_samples; g++)
            {
                for (uint32_t p = 0; p < pack.s.ploidy; p++)
                {
                    start =  (g*pack.s.ploidy+p);
                    pt[start] = lut[resAll[start]][resAll[start+no_haplotypes]][(perm[sampleIDs[g]*pack.s.ploidy + p])%8];
                    //memcpy(pt + (start), lut[resAll[start]][resAll[start+no_haplotypes]] +(perm[sampleIDs[g]*pack.s.ploidy + p])%8, 1);
                }
            }
            
            str.l = str.l + (smpl.no_samples*pack.s.ploidy);
            str.s[str.l] = 0;
            
            bcf_update_info_int32(hdr, record, "_row", NULL, 0);
            // AC/AN count
            if(out_AC_AN)
            {
                if(setACAN(hdr, record, str))
                {
                    if(out_genotypes)
                        bcf_update_genotypes_fast(str, record);
                    
                    bcf_write1(out, hdr, record);
                    written_records++;
                    
                }
            }
            else
            {
                if(out_genotypes)
                    bcf_update_genotypes_fast(str, record);
                
                bcf_write1(out, hdr, record);
                written_records++;
            }
            
            var_idx++;
            
            
            prev_block_id = block_id;
        }
    }
    hts_close(out);
    free(tmpia);
    if(str.s)
    {
        free(str.s);
        str.s = nullptr;
    }
    bcf_destroy(record);

    if(resUnique)
        delete [] resUnique;
    if(resAll)
        delete [] resAll;
    if(perm)
        delete [] perm;
    return 0;
}

int Decompressor::decompressRangeSample(const string & range)
{
    initialLut();
    
    uchar_t *decomp_data = nullptr;
    uchar_t * decomp_data_perm = nullptr;
    uint32_t * perm = nullptr;
    uint32_t * rev_perm = nullptr;
    uint32_t pos;
    
    uint32_t no_haplotypes = smpl.no_samples*pack.s.ploidy;
    perm = new uint32_t[pack.s.n_samples * pack.s.ploidy];
    rev_perm = new uint32_t[pack.s.n_samples * pack.s.ploidy];
    
    bcf1_t * record = bcf_init();
    uint32_t gg, vec1_start, vec2_start;
//    uint32_t  end;
    uint32_t i = 0;
    
    buff_bm.setBitMemory(&pack.bm);
    
    const char *key = "GT";
    khint_t k;
    vdict_t *d;
    int fmt_id;
    done_unique.clear();
    uint32_t block_id, prev_block_id = 0xFFFF;
    
    uint32_t written_records = 0;
    kstring_t str = {0,0,0};
    
    char * tmp_vec = new char[pack.s.vec_len * 8 + 8];

    long long * tmp_vec_ll = (long long *)tmp_vec;
    long long * lut_ll = nullptr;
    
    
    if(range == "")
    {
        d = (vdict_t*)hdr->dict[BCF_DT_ID];
        k = kh_get(vdict, d, key);
        fmt_id = (k == kh_end(d)? -1 : kh_val(d, k).id);
        
        if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,fmt_id) )
        {
            if(perm)
                delete [] perm;
            if(rev_perm)
                delete [] rev_perm;
            bcf_destroy(record);
            if ( !no_haplotypes ) return 0;
            return -1;  // the key not present in the header
        }
        
        decomp_data = new uchar_t[pack.s.vec_len*2];
        
        decomp_data_perm = new uchar_t[pack.s.vec_len*2];

        bcf_enc_int1(&str, fmt_id);
        
        bcf_enc_size(&str, pack.s.ploidy, BCF_BT_INT8);
        
        char *tmp;
        str.m = str.l + no_haplotypes + 1;
        kroundup32(str.m);
        if ((tmp = (char*)realloc(str.s, str.m)))
            str.s = tmp;
        else
            exit(8);
        
        str.l = 3;
        if (no_haplotypes == 0) bcf_enc_size(&str, 0, BCF_BT_NULL);
        
       uint32_t  end;
        if(no_haplotypes & 7)//%8)
        {
            end = pack.s.vec_len - 1;
            gg = no_haplotypes & 7;//%8;
        }
        else
        {
            end = pack.s.vec_len;
            gg = 0;
        }
        
        
        
        while (bcf_read1(bcf, hdr, record)>= 0 && written_records < records_to_process)
        {
            pos = 0;
            str.l = 3;
            
            bcf_unpack(record, BCF_UN_ALL);
            record->n_sample = bcf_hdr_nsamples(hdr);
            
            // Permutations
            block_id = i/pack.s.max_no_vec_in_block; //pair of vectors always in the same block (so it is enough to find first vector block
            if(block_id != prev_block_id)
            {
                // Get approproate permutations
                pack.getPermArray(block_id, perm);
                reverse_perm(perm, rev_perm, pack.s.n_samples*pack.s.ploidy);
            }
            
            fill_n(decomp_data, pack.s.vec_len*2, 0);
            
            vec1_start = 0;
            vec2_start = pack.s.vec_len;
            decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
            decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
            
            // Permutations
            decode_perm_rev(pack.s.n_samples*pack.s.ploidy, vec2_start, rev_perm, decomp_data_perm, decomp_data);
            
            prev_block_id = block_id;
            
            char *pt = str.s + str.l;
         /*   for (uint32_t g = 0; g < smpl.no_samples; ++g)
            {
                for(int p = 0; p < (int) pack.s.ploidy; p++)
                {
                    vec1_start = (sampleIDs[g]*pack.s.ploidy + p)/8; //2 vectors per variant
                    memcpy(pt + (g * pack.s.ploidy + p), lut[decomp_data[vec1_start]][decomp_data[vec1_start + pack.s.vec_len]] + (sampleIDs[g]*pack.s.ploidy + p)%8, 1);
                }
            }
            
            str.l = str.l + smpl.no_samples * pack.s.ploidy;
            str.s[str.l] = 0;*/
            
            
            /////////////////

          //  uint32_t start;
            for (vec1_start = 0; vec1_start < end; ++vec1_start)
            {
                //memcpy(tmp_vec + (vec1_start << 3), lut[decomp_data[vec1_start]][decomp_data[vec2_start++]], 8);
                lut_ll = (long long *)(lut[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                *(tmp_vec_ll + vec1_start) = *lut_ll;//*(long long *)lut[decomp_data[vec1_start]][decomp_data[vec2_start++]];
            }
            
            // str.l = str.l + (end << 3);
            if(gg)
            {
                memcpy(tmp_vec + (end << 3), lut[decomp_data[vec1_start]][decomp_data[vec2_start]], gg);
                //      str.l = str.l + g;
            }
            
            
            uint32_t which;
             if(pack.s.ploidy == 2)
             {
                for (uint32_t g = 0; g < smpl.no_samples ; ++g)
                {
                    which = sampleIDs[g]*pack.s.ploidy;
                    pt[g * pack.s.ploidy] = tmp_vec[which];
                    pt[g * pack.s.ploidy + 1] = tmp_vec[which + 1];
                }
             }
             else if (pack.s.ploidy == 1)
             {
                for (uint32_t g = 0; g < smpl.no_samples ; ++g)
                {
                    which = sampleIDs[g]*pack.s.ploidy;
                    pt[g * pack.s.ploidy] = tmp_vec[which];
                }
             }
             else
             {
                for (uint32_t g = 0; g < smpl.no_samples ; ++g)
                {
                    which = sampleIDs[g]*pack.s.ploidy;
                    for(int p = 0; p < (int) pack.s.ploidy; p++)
                        pt[g * pack.s.ploidy + p] = tmp_vec[which + p];
                }
             }
            
            
            
             /*
            for (uint32_t g = 0; g < smpl.no_samples ; ++g)
            {
               // for(int p = 0; p < (int) pack.s.ploidy; p++)
                {
                    which = sampleIDs[g]*pack.s.ploidy;
                    memcpy(pt + (g * pack.s.ploidy), tmp_vec + which, pack.s.ploidy);
                }
            }*/
            
            
            
            str.l = str.l + smpl.no_samples*pack.s.ploidy;
            str.s[str.l] = 0;
            
            bcf_update_info_int32(hdr, record, "_row", NULL, 0);
       
            // AC/AN count
            if(out_AC_AN)
            {
                if(setACAN(hdr, record, str))
                {
                    if(out_genotypes)
                        bcf_update_genotypes_fast(str, record);
                    
                    bcf_write1(out, hdr, record);
                    written_records++;
                    
                }
            }
            else
            {
                if(out_genotypes)
                    bcf_update_genotypes_fast(str, record);
                
                bcf_write1(out, hdr, record);
                written_records++;
            }
            
        }
        
        for (auto & it:done_unique)
        {
            delete [] it.second;
        }
        done_unique.clear();
        
        delete [] tmp_vec;

    }
    else
    {
        d = (vdict_t*)hdr->dict[BCF_DT_ID];
        k = kh_get(vdict, d, key);
        fmt_id = (k == kh_end(d)? -1 : kh_val(d, k).id);
        
        if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,fmt_id) )
        {
            if(perm)
                delete [] perm;
            if(rev_perm)
                delete [] rev_perm;
            bcf_destroy(record);
            
            if ( !no_haplotypes ) return 0;
            return -1;  // the key not present in the header
        }
        
        decomp_data = new uchar_t[pack.s.vec_len*2];
        decomp_data_perm = new uchar_t[pack.s.vec_len*2];

        
        bcf_enc_int1(&str, fmt_id);
        bcf_enc_size(&str, pack.s.ploidy, BCF_BT_INT8);

        char *tmp;
        str.m = str.l + no_haplotypes + 1;
        kroundup32(str.m);
        if ((tmp = (char*)realloc(str.s, str.m)))
            str.s = tmp;
        else
            exit(8);
        str.l = 3;
        if (no_haplotypes == 0) bcf_enc_size(&str, 0, BCF_BT_NULL);
        
        uint32_t  end;
        if(no_haplotypes & 7)//%8)
        {
            end = pack.s.vec_len - 1;
            gg = no_haplotypes & 7;//%8;
        }
        else
        {
            end = pack.s.vec_len;
            gg = 0;
        }

        
        
        
        bcf_info_t * a;
        hts_itr_t * itr = bcf_itr_querys(bcf_idx, hdr, range.c_str());
        
        char * tmp_vec = new char[pack.s.vec_len * 8 + 8];
        
    
        
        while(bcf_itr_next(bcf, itr, record) != -1 && written_records < records_to_process)
        {
            if(!i)
            {
                a = bcf_get_info(hdr, record, "_row"); //bcf_get_info_id
                i = (a->v1.i*2);
            }
            
            pos = 0;
            str.l = 3;
            
            bcf_unpack(record, BCF_UN_ALL);
            record->n_sample = bcf_hdr_nsamples(hdr);
            
            // Permutations
            block_id = i/pack.s.max_no_vec_in_block; //pair of vectors always in the same block (so it is enough to find first vector block
            if(block_id != prev_block_id)
            {
                // Get approproate permutations
                pack.getPermArray(block_id, perm);
                reverse_perm(perm, rev_perm, pack.s.n_samples*pack.s.ploidy);
            }
            fill_n(decomp_data, pack.s.vec_len*2, 0);
            
            vec1_start = 0;
            vec2_start = pack.s.vec_len;
            decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
            decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
            
            // Permutations
          //  decode_perm(pack.s.n_samples * pack.s.ploidy, vec2_start, perm, decomp_data_perm, decomp_data);
           
            decode_perm_rev(pack.s.n_samples*pack.s.ploidy, vec2_start, rev_perm, decomp_data_perm, decomp_data);
            prev_block_id = block_id;
            
            char *pt = str.s + str.l;
         /*
            for (uint32_t g = 0; g < smpl.no_samples ; ++g)
            {
                for(int p = 0; p < (int) pack.s.ploidy; p++)
                {
                    vec1_start = (sampleIDs[g]*pack.s.ploidy + p)/8; //2 vectors per variant
                    memcpy(pt + (g * pack.s.ploidy + p), lut[decomp_data[vec1_start]][decomp_data[vec1_start + pack.s.vec_len]]+(sampleIDs[g] * pack.s.ploidy + p)%8, 1);
                }
            }
            
            str.l = str.l + smpl.no_samples*pack.s.ploidy;
            str.s[str.l] = 0;
          // cout << i << " " << str.l << " " <<  (int) str.s[0] << (int) str.s[1] << (int) str.s[2] << (int) str.s[3] << (int) str.s[4] << (int) str.s[2002] << endl;
          */
            
            
            
             /////////////////

            for (vec1_start = 0; vec1_start < end; ++vec1_start)
            {
               // memcpy(tmp_vec + (vec1_start << 3), lut[decomp_data[vec1_start]][decomp_data[vec2_start++]], 8);
                lut_ll = (long long *)(lut[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                *(tmp_vec_ll + vec1_start) = *lut_ll;//*(long long *)lut[decomp_data[vec1_start]][decomp_data[vec2_start++]];

            }
         
           // str.l = str.l + (end << 3);
            if(gg)
            {
                memcpy(tmp_vec + (end << 3), lut[decomp_data[vec1_start]][decomp_data[vec2_start]], gg);
          //      str.l = str.l + g;
            }
  
            
            uint32_t which;
            if(pack.s.ploidy == 2)
            {
                for (uint32_t g = 0; g < smpl.no_samples ; ++g)
                {
                    which = sampleIDs[g]*pack.s.ploidy;
                    pt[g * pack.s.ploidy] = tmp_vec[which];
                    pt[g * pack.s.ploidy + 1] = tmp_vec[which + 1];
                }
            }
            else if (pack.s.ploidy == 1)
            {
                for (uint32_t g = 0; g < smpl.no_samples ; ++g)
                {
                    which = sampleIDs[g]*pack.s.ploidy;
                    pt[g * pack.s.ploidy] = tmp_vec[which];
                }
            }
            else
            {
                for (uint32_t g = 0; g < smpl.no_samples ; ++g)
                {
                    which = sampleIDs[g]*pack.s.ploidy;
                    for(int p = 0; p < (int) pack.s.ploidy; p++)
                        pt[g * pack.s.ploidy + p] = tmp_vec[which + p];
                }
            }

            
            
            str.l = str.l + smpl.no_samples*pack.s.ploidy;
            str.s[str.l] = 0;
            
            
            bcf_update_info_int32(hdr, record, "_row", NULL, 0);
       
            // AC/AN count
            if(out_AC_AN)
            {
                if(setACAN(hdr, record, str))
                {
                    if(out_genotypes)
                        bcf_update_genotypes_fast(str, record);
                    
                    bcf_write1(out, hdr, record);
                    written_records++;
                    
                }
            }
            else
            {
                if(out_genotypes)
                    bcf_update_genotypes_fast(str, record);
                
                bcf_write1(out, hdr, record);
                written_records++;
            }
            
        }
        
        bcf_itr_destroy(itr);
        
        for (auto & it:done_unique)
        {
            delete [] it.second;
        }
        done_unique.clear();
        
        delete [] tmp_vec;
    }
    
    hts_close(out);
    
    if(decomp_data)
        delete [] decomp_data;
    if(decomp_data_perm)
        delete [] decomp_data_perm;
    if(perm)
        delete perm;
    if(rev_perm)
       delete[] rev_perm;
    
    if(str.s)
    {
        free(str.s);
        str.s = nullptr;
    }
    bcf_destroy(record);
    return 0;
}


bool Decompressor::setACAN(bcf_hdr_t * hdr, bcf1_t * record, const kstring_t & str)
{
    int32 an = 0, ac = 0;

    for(int32_t l = 3; l < (int) str.l; l++)
    {
        an++;
        if(str.s[l] == 0x05)
            ac++;
    }

    if(ac >= (int) minAC && ac <= (int) maxAC)
    {
        
        bcf_update_info_int32(hdr, record, "AN", &an, 1);
        bcf_update_info_int32(hdr, record, "AC", &ac, 1);
        return true;
    }

    return false;
}

/************************/
// full_decode: if true, always decode full; otherwise decode only unique vectors
uchar_t Decompressor::get_vec_bytes(uint64_t vec_id, std::vector< std::pair<uint32_t, uint32_t> > & whichByte_whereInRes, uchar_t * resUnique, bool & is_uniqe_id, uint64_t & curr_zeros, uint64_t & curr_copy, uchar_t * resAll, uint32_t unique_pos_first_in_block, uint32_t first_vec_in_block, bool full_decode) //uchar_t Decompressor::get_vec_bytes(uint64_t vec_id, std::vector< std::pair<uint32_t, uchar_t> > & whichByte_whereInRes, uchar_t * resUnique, bool & is_uniqe_id, uint64_t & curr_zeros, uint64_t & curr_copy, uchar_t * resAll)
{
    
   //if(vec_id == 322560)
   //     cout << "here" << endl;

    
    uint32_t next_haplotype = 0; // Next byte to decode
    
    uint32_t last_byte = whichByte_whereInRes.back().first; // Last byte to decode
    
    uchar_t byte;
    
    uint32_t tmp = 0;
    
    uint8_t parity = vec_id&1;
    uint64_t vector = vec_id >> 1, curr_non_copy_vec_id;
    uint64_t no_haplotypes = whichByte_whereInRes.size() - 1;
    uint64_t vec_start = parity*no_haplotypes;//(vec_id - first_vec_in_block) * no_haplotypes;
    
    if(parity)
    {
        if(copy_bit_vector[parity][vector]) // If vector is a copy od other vector (certainly it is placed within the same block)
        {
            
            if(!full_decode)
            {
                is_uniqe_id = false;
                curr_copy++;
                return 0;
            }
            
            unsigned long long bit_pos = (curr_copy)*pack.used_bits_cp;
            
            pack.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
            pack.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos&7); // %8
            pack.bm_comp_copy_orgl_id.GetBits(tmp, pack.used_bits_cp);
            
            // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
            curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy;
            
            curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;
            is_uniqe_id = false;
            curr_copy++;
            
            
           
          //  for(uint32_t idx = 0; idx < no_haplotypes; idx++)
            {
               // resAll[whichByte_whereInRes[idx].second * pack.no_vec + vec_id] = resUnique[pack.no_non_copy*whichByte_whereInRes[idx].second  + curr_non_copy_vec_id];
             //   resAll[vec_id * pack.s.vec_len + whichByte_whereInRes[idx].second] = resUnique[curr_non_copy_vec_id*pack.s.vec_len + whichByte_whereInRes[idx].second];
                memcpy(resAll + vec_start, resUnique + (curr_non_copy_vec_id-unique_pos_first_in_block)*no_haplotypes, no_haplotypes) ;
                
            }
            return 0;
        }
        else if(zeros_only_bit_vector[1][vector])
        {
            is_uniqe_id = false;
            curr_zeros++;
            
            if(!full_decode)
                return 0;
            
          //  for(uint32_t idx = 0; idx < no_haplotypes; idx++)
            {
//                resAll[idx * pack.no_vec + vec_id] = 0;  //zero everywhere
            //    resAll[vec_id * pack.s.vec_len + idx] = 0;
                
                fill_n(resAll + vec_start, no_haplotypes, 0);

                
            }
            
            return 0;
        }
        else
        {
            curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy;
            is_uniqe_id = true;
        }
    }
    else if(!zeros_only_bit_vector[0][vector])
    {
        is_uniqe_id = false;
        curr_zeros++;
        
        if(!full_decode)
            return 0;
        
      /*  for(uint32_t idx = 0; idx < no_haplotypes; idx++)
        {
           // resAll[idx * pack.no_vec + vec_id] = 0;
            resAll[vec_id * pack.s.vec_len + vec_id] = 0;
        }*/
        fill_n(resAll + vec_start, no_haplotypes, 0);
        
        return 0;
    }
    else
    {
        if(copy_bit_vector[parity][vector])
        {
            if(!full_decode)
            {
                is_uniqe_id = false;
                curr_copy++;
                return 0;
            }
            
            unsigned long long bit_pos = (curr_copy)*pack.used_bits_cp;
            
            pack.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
            pack.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos&7); // %8
            pack.bm_comp_copy_orgl_id.GetBits(tmp, pack.used_bits_cp);
            
            // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
            curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy;
            curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;
            is_uniqe_id = false;
            curr_copy++;
            
            
            if(!full_decode)
                return 0;
        //    for(uint32_t idx = 0; idx < no_haplotypes; idx++)
            {
             //   resAll[whichByte_whereInRes[idx].second * pack.no_vec + vec_id] = resUnique[pack.no_non_copy*whichByte_whereInRes[idx].second  + curr_non_copy_vec_id];
                
          //      resAll[vec_id * pack.s.vec_len + whichByte_whereInRes[idx].second] = resUnique[curr_non_copy_vec_id*pack.s.vec_len + whichByte_whereInRes[idx].second];
                
                memcpy(resAll + vec_start, resUnique + (curr_non_copy_vec_id-unique_pos_first_in_block)*no_haplotypes, no_haplotypes);
                
            }
            return 0;
        }
        else
        {
            curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy;
            is_uniqe_id = true;
        }
    }
    
    uint32_t curr_pos;
    tmp = 0;
    unsigned long long full_pos = curr_non_copy_vec_id/FULL_POS_STEP  * (sizeof(uint32_t) + pack.used_bits_noncp*(FULL_POS_STEP-1)/BITS_IN_BYTE);
    pack.bm_comp_pos.SetPos(full_pos);
    pack.bm_comp_pos.GetWord(curr_pos);
    uint32_t j = ((curr_non_copy_vec_id-1)%FULL_POS_STEP)/(sizeof(uint32_t)*BITS_IN_BYTE), end = curr_non_copy_vec_id%FULL_POS_STEP;
    if(j)
    {
        pack.bm_comp_pos.SetPos(full_pos + sizeof(uint32_t) + j*pack.used_bits_noncp*sizeof(uint32_t));
        j = j*(sizeof(uint32_t)*BITS_IN_BYTE);
        
        if(end > j + 1)
            pack.bm_comp_pos.GetBitsAndDiscard((int32_t)pack.used_bits_noncp * (end - 1 - j));
        if(end  > j)
            pack.bm_comp_pos.GetBits(tmp, (int32_t)pack.used_bits_noncp);
    }
    else
    {
        if(end > 1)
            pack.bm_comp_pos.GetBitsAndDiscard((int32_t)pack.used_bits_noncp*(end - 1));
        if(end)
            pack.bm_comp_pos.GetBits(tmp, (int32_t)pack.used_bits_noncp);
        
    }
    curr_pos += tmp;
    
    uint32_t ones_group;
    
    uint32_t decoded_bytes = 0;
    uint32_t zero_run_len,ones_run_len;
    uint32_t best_pos = 0;
    uint32_t best_match_len = 0;
    int32_t  flag;
    uint32 litRun;
    
    buff_bm.SetPos(curr_pos);
    
    ones_group = buff_bm.decodeFastLut(&pack.huf_group_type);
    
    CHuffman * h_lit = &pack.huf_literals[ones_group];
    
    while(decoded_bytes <= last_byte)
    {
        flag = buff_bm.decodeFastLut(&pack.huf_flags);
        
        switch(flag)
        {
            case 0: //literal
            {
                byte = buff_bm.decodeFastLut(h_lit);
                while(decoded_bytes == whichByte_whereInRes[next_haplotype].first) //curr byte == next byte to decode
                {
                  //  resAll[whichByte_whereInRes[next_haplotype].second * pack.no_vec + vec_id] = byte;
                    resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = byte;
                    next_haplotype++;
//                    if(next_haplotype == no_haplotypes)
//                        return 0; //all haplotypes decoded
                }
                
                decoded_bytes++;
                break;
            }
            case 1: //match
            {
                // Difference between current and match is stored
                tmp = buff_bm.decodeFastLut(&pack.huf_match_diff_MSB);
                best_pos = tmp << (pack.s.bit_size_id_match_pos_diff - MATCH_BITS_HUF);
                tmp = buff_bm.getBits(pack.s.bit_size_id_match_pos_diff - MATCH_BITS_HUF);
                best_pos = best_pos | tmp;
                
                best_pos += 1;
                best_pos = curr_non_copy_vec_id - best_pos;
                
                //  prev_vec_match = best_pos;
                best_match_len = buff_bm.decodeFast(&pack.huf_match_lens[ones_group]);
                
                while(decoded_bytes + best_match_len  >  whichByte_whereInRes[next_haplotype].first)
                {
                  //  resAll[whichByte_whereInRes[next_haplotype].second * pack.no_vec + vec_id] = resUnique[pack.no_non_copy*whichByte_whereInRes[next_haplotype].second  + best_pos];
                    resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = resUnique[(best_pos-unique_pos_first_in_block)*no_haplotypes + whichByte_whereInRes[next_haplotype].second];
                    
                    next_haplotype++;
//                    if(next_haplotype == no_haplotypes)
//                        return 0; //all haplotypes decoded
                }
                decoded_bytes += best_match_len;
                break;
            }
            case 2: //match same
            {
                best_match_len = buff_bm.decodeFast(&pack.huf_match_lens[ones_group]);
               
                while(decoded_bytes + best_match_len  >  whichByte_whereInRes[next_haplotype].first)
                {
                  //  resAll[whichByte_whereInRes[next_haplotype].second * pack.no_vec + vec_id] = resUnique[pack.no_non_copy*whichByte_whereInRes[next_haplotype].second  + best_pos];
                     resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = resUnique[(best_pos-unique_pos_first_in_block)*no_haplotypes + whichByte_whereInRes[next_haplotype].second];
                    
                    next_haplotype++;
//                    if(next_haplotype == no_haplotypes)
//                        return 0; //all haplotypes decoded
                }
                
                decoded_bytes += best_match_len;
                break;
            }
            case 3:
            {
                zero_run_len = buff_bm.decodeFast(&(pack.huf_zeros_runs[ones_group]));
                
                while(decoded_bytes + zero_run_len  >  whichByte_whereInRes[next_haplotype].first)
                {
                    //resAll[whichByte_whereInRes[next_haplotype].second * pack.no_vec + vec_id] = 0;
                    
                    resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = 0;
                    
                    next_haplotype++;
//                    if(next_haplotype == no_haplotypes)
//                        return 0; //all haplotypes decoded
                }
                
                decoded_bytes += zero_run_len;
                break;
            }
            case 4:
            {
                ones_run_len = buff_bm.decodeFast(&(pack.huf_ones_runs[ones_group]));
                
                while(decoded_bytes + ones_run_len  >  whichByte_whereInRes[next_haplotype].first)
                {
                   // resAll[whichByte_whereInRes[next_haplotype].second * pack.no_vec + vec_id] = 0xFF;
                    resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = 0xFF;
                    
                    next_haplotype++;
//                    if(next_haplotype == no_haplotypes)
//                        return 0; //all haplotypes decoded
                }
                
                decoded_bytes += ones_run_len;
                break;
            }
            default: //run of literals (2 (MIN_LITERAL_RUN) - MAX_LITERAL_RUN)
            {
                flag = flag - 3; // Shift by 3, to be able to difference between 2,3 and 4 flags nad literal run of 2,3 or 4
                
                if(flag + decoded_bytes <= whichByte_whereInRes[next_haplotype].first)
                {
                    litRun = buff_bm.getBits(pack.used_bits_litRunSize[ones_group]);
                    if(!litRun) //litRun == 0 means used_bits_litRunSize bits were not enough to store size
                        litRun = buff_bm.getBits(pack.max_used_bits_litRunSize[ones_group]);
                    
                    litRun = litRun + pack.minLitRunSize[ones_group][flag] - 1; //1 is added to litRun, so it is never == 0
                    
                    //skip run of literals
                    buff_bm.getBitsAndDiscard(litRun);
                    decoded_bytes += flag;
                }
                else
                {
                    litRun = buff_bm.getBits(pack.used_bits_litRunSize[ones_group]);
                    if(!litRun) //litRun == 0 means used_bits_litRunSize bits were not enough to store size
                        buff_bm.getBitsAndDiscard(pack.max_used_bits_litRunSize[ones_group]);

                    for(int i = 0; i < flag; i++)
                    {
                        byte = buff_bm.decodeFastLut(h_lit);
                        while(decoded_bytes  ==  whichByte_whereInRes[next_haplotype].first)
                        {
                            //resAll[whichByte_whereInRes[next_haplotype].second * pack.no_vec + vec_id] = byte;
                            
                            resAll[vec_start + whichByte_whereInRes[next_haplotype].second] =  byte;
                            next_haplotype++;
//                            if(next_haplotype == no_haplotypes)
//                                return 0; //all haplotypes decoded
                        }
                        
                        decoded_bytes++;
                    }
                }
            }
        }
        
        if(next_haplotype == no_haplotypes)
            return 0; //all haplotypes decoded
        
    }
    return 0;
}

bool Decompressor::loadPack()
{
    bool b = pack.loadPack(arch_name);
    if(b)
    {
        zeros_only_vector = new uchar_t[pack.s.vec_len]();
        ones_only_vector = new uchar_t[pack.s.vec_len];
        fill_n(ones_only_vector, pack.s.vec_len, 0xFF);
    }
    return b;
}

int Decompressor::initOut()
{
    char write_mode[5] = "wb-";
    if(out_type == BCF)
    {
        write_mode[3] = compression_level;
        write_mode[4] = '\0';
    }
    if(out_name != "")
    {
        char *gz_fname = (char*) malloc(strlen(out_name.c_str())+5);
        if(out_type == VCF)
        {
            snprintf(gz_fname,strlen(out_name.c_str())+5,"%s.vcf",out_name.c_str());
            out   = hts_open(gz_fname, "w");
        }
        else
        {
            snprintf(gz_fname,strlen(out_name.c_str())+5,"%s.bcf",out_name.c_str());
            out   = hts_open(gz_fname, write_mode);
        }
        free(gz_fname);
    }
    else
    {
        if(out_type == VCF)
        {
            out = hts_open("-", "w");
        }
        else
        {
            out = hts_open("-", write_mode);
        }
        
    }
    
    string filename(arch_name);
    filename = filename + ".ind";
    
    if(samples_to_decompress == "")
        smpl.setAllSamples(hdr, filename, out_genotypes);
    else
    {
        smpl.loadSamples(filename);
        sampleIDs = smpl.setSamples(hdr, samples_to_decompress.c_str(), out_genotypes);
    }
    
    bcf_hdr_add_sample(hdr, NULL);
    bcf_hdr_write(out, hdr);
    
    if(out_AC_AN)
    {
        bcf_hdr_append(hdr,"##INFO=<ID=AC,Number=A,Type=String,Description=\"Count of alternate alleles\">");
        bcf_hdr_append(hdr,"##INFO=<ID=AN,Number=A,Type=String,Description=\"Count of total alleles\">");
        bcf_hdr_sync(hdr);
        
        //more strict of minAC/minAF is taken into account
        if(minAF && minAF > (double)minAC/(smpl.no_samples*pack.s.ploidy))
            minAC = ceil(minAF * (smpl.no_samples*pack.s.ploidy));
        //more strict of minAC/minAF is taken into account
        if(maxAF < 1 && maxAF < (double)maxAC/(smpl.no_samples*pack.s.ploidy))
            maxAC = floor(maxAF * (smpl.no_samples*pack.s.ploidy));
        //minAC and maxAC are used
    }
    return 0;
}

int Decompressor::loadBCF()
{
    char *fname = (char*) malloc(arch_name.length()+5);
    char *fname2 = (char*) malloc(arch_name.length()+10);
    // BCF load
    snprintf(fname,arch_name.length() + 5,"%s.bcf",arch_name.c_str());
    bcf = hts_open(fname, "rb");
    snprintf(fname2,arch_name.length() + 9,"%s.bcf.csi",arch_name.c_str());
    bcf_idx = hts_idx_load2(fname, fname2);
    hdr = bcf_hdr_read(bcf);
    
    free(fname);
    free(fname2);
    return 0;
}

void Decompressor::setMemoryUsage()
{
    // Set memory usage
    if(MB_memory)
    {
        if(max_MB_memory)
            max_stored_unique = (max_MB_memory*1000000) / pack.s.vec_len;
        else
            max_stored_unique = pack.no_vec;
        
        if (PART_SIZE*2 < max_stored_unique)
            max_stored_unique = PART_SIZE*2;
    }
    else
        max_stored_unique = 0;
}

void Decompressor::initialLut()
{
    uchar mask;
    
    for(int c = 0; c < 8; c++)
    {  mask = 0x80 >> c;
        for(int i = 0; i < 256; i++)
        {
            if(i & mask)
            {
                for(int j = 0; j < 256; j++)
                {
                    if(j & mask)  //11
                    {lut[i][j][c] = bcf_gt_phased(2);}
                    else            //10
                    {lut[i][j][c] = bcf_gt_missing;}
                }
                
            }
            else
            {
                for(int j = 0; j < 256; j++)
                {
                    if(j & mask)    //01
                    {lut[i][j][c] = bcf_gt_phased(1);}
                    else                //00
                    {lut[i][j][c] = bcf_gt_phased(0);}
                }
            }
        }
    }
	
	// Permutation LUT
	for(int i = 0; i < 8; ++i)
		perm_lut[i] = 1 << (7 - i);	
}

void Decompressor::decomp_vec_rrr_range(uint64_t vec_id, uint64_t offset, uint64_t length, uint32_t & pos, uchar_t *decomp_data, uint64_t start_id, bool is_unique_id)
{
    
    uint64_t id, curr_non_copy_vec_id, toDelete;
    uint32_t curr_pos, decoded_bytes;
    uint32_t zero_run_len, ones_run_len, rest=0;
    uint32_t best_pos = 0;
    uint32_t best_match_len = 0;
    uint64_t prev_vec_match = 0;
    int32_t  flag;
    uint32_t bit0;
    int32 bits_read;
    uint32 tmp = 0;
    uint16_t left;
    
    if(is_unique_id)
    {
        curr_non_copy_vec_id = vec_id;
        
        got_it = done_unique.find (curr_non_copy_vec_id);
        
        if ( got_it != done_unique.end())
        {
            memcpy(decomp_data+pos, got_it->second + offset, length);
            pos += length;
            return;
        }
        
    }
    else
    {
        id = vec_id >> 1; //vec_id/2
        uchar parity = vec_id&1;
        if((parity && pack.rrr_zeros_only_bit_vector[1][id]) || (!(parity) && !pack.rrr_zeros_only_bit_vector[0][id]))
        {
            memcpy(decomp_data+pos, zeros_only_vector, pack.s.vec_len);
            pos += pack.s.vec_len;
            return;
        }
        
        if(pack.rrr_copy_bit_vector[parity][id])
        {
            unsigned long long bit_pos = (pack.rrr_rank_copy_bit_vector[0](id + ((parity))) + pack.rrr_rank_copy_bit_vector[1](id))*pack.used_bits_cp;
            
            pack.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);  // /8
            pack.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos&7);  // %8
            pack.bm_comp_copy_orgl_id.GetBits(tmp, pack.used_bits_cp);
            
            // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
            curr_non_copy_vec_id = vec_id - pack.rrr_rank_zeros_only_bit_vector_0(id+((parity))) - pack.rrr_rank_zeros_only_bit_vector_1(id) - \
            pack.rrr_rank_copy_bit_vector[0](id+((parity))) - pack.rrr_rank_copy_bit_vector[1](id);
            
            curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;
            
            got_it = done_unique.find (curr_non_copy_vec_id);
            
            if ( got_it != done_unique.end())
            {
                memcpy(decomp_data+pos, got_it->second + offset, length);
                pos += length;
                return;
            }
        }
        else //unique and not a copy - no need to check if it was previously decompressed (got_it)
        {
            curr_non_copy_vec_id = vec_id - pack.rrr_rank_zeros_only_bit_vector_0(id+((parity))) - pack.rrr_rank_zeros_only_bit_vector_1(id) - \
            pack.rrr_rank_copy_bit_vector[0](id+((parity))) - pack.rrr_rank_copy_bit_vector[1](id);
        }
    }

    tmp = 0;
    unsigned long long full_pos = curr_non_copy_vec_id/FULL_POS_STEP  * (sizeof(uint32_t) + pack.used_bits_noncp*(FULL_POS_STEP-1)/BITS_IN_BYTE);
    pack.bm_comp_pos.SetPos(full_pos);
    pack.bm_comp_pos.GetWord(curr_pos);
    uint32_t  j = ((curr_non_copy_vec_id-1)%FULL_POS_STEP)/(sizeof(uint32_t)*BITS_IN_BYTE), end = curr_non_copy_vec_id%FULL_POS_STEP;
    if(j)
    {
        pack.bm_comp_pos.SetPos(full_pos + sizeof(uint32_t) + j*pack.used_bits_noncp*sizeof(uint32_t));
        j = j*(sizeof(uint32_t)*BITS_IN_BYTE);
        
        if(end > j + 1)
            pack.bm_comp_pos.GetBitsAndDiscard((int32_t)pack.used_bits_noncp * (end - 1 - j));
        if(end  > j)
            pack.bm_comp_pos.GetBits(tmp, (int32_t)pack.used_bits_noncp);
    }
    else
    {
        if(end > 1)
            pack.bm_comp_pos.GetBitsAndDiscard((int32_t)pack.used_bits_noncp*(end - 1));
        if(end)
            pack.bm_comp_pos.GetBits(tmp, (int32_t)pack.used_bits_noncp);
    }
    curr_pos += tmp;
    
    uint32 litRun;
    buff_bm.SetPos(curr_pos);
    uint32_t ones_group;
    ones_group = buff_bm.decodeFastLut(&pack.huf_group_type);
    CHuffman * h_lit = &pack.huf_literals[ones_group];
    decoded_bytes = 0;
    while(decoded_bytes < offset)
    {
        flag = buff_bm.decodeFastLut(&pack.huf_flags);
        switch(flag)
        {
            case 0:  //literal
            {
                buff_bm.decodeFastLut(h_lit);
                decoded_bytes++;
                break;
            }
            case 1: //match
            {
                // Difference between current and match id
                tmp = buff_bm.decodeFastLut(&pack.huf_match_diff_MSB);
                best_pos = tmp << (pack.s.bit_size_id_match_pos_diff - MATCH_BITS_HUF);
                tmp = buff_bm.getBits(pack.s.bit_size_id_match_pos_diff - MATCH_BITS_HUF);
                best_pos = best_pos | tmp;
                best_pos += 1; //shift (to not waste 1 value)
                
                best_pos = curr_non_copy_vec_id - best_pos;
                prev_vec_match = best_pos;
                
                best_match_len = buff_bm.decodeFast(&pack.huf_match_lens[ones_group]);
                if(decoded_bytes + best_match_len  <=  offset)
                {
                    decoded_bytes += best_match_len;
                }
                else if(decoded_bytes + best_match_len <= offset + length)
                {
                    buff_bm.getBuffer(bit0, left);
                    curr_pos = pack.bm.GetPos() - 1;
                    bits_read = 8 - pack.bm.GetWordPos();
                    nesting++;
                    decomp_vec_rrr_range(best_pos, offset, decoded_bytes + best_match_len - offset, pos, decomp_data, start_id, true);
                    nesting--;
                    pack.bm.SetPos(curr_pos);
                    pack.bm.GetBitsAndDiscard(bits_read);
                    
                    buff_bm.setBuffer(bit0, left);
                    rest = decoded_bytes + best_match_len  - offset;
                    decoded_bytes = offset;
                    
                }
                else //if(decoded_bytes + best_match_len > offset + length)
                {
                    buff_bm.getBuffer(bit0, left);
                    
                    curr_pos = pack.bm.GetPos() - 1;
                    bits_read = 8 - pack.bm.GetWordPos();
                    nesting++;
                    decomp_vec_rrr_range(best_pos, offset, length, pos, decomp_data, start_id, true);
                    nesting--;
                    pack.bm.SetPos(curr_pos);
                    pack.bm.GetBitsAndDiscard(bits_read);
                    buff_bm.setBuffer(bit0, left);
                    decoded_bytes = offset;
                    rest = length;
                }
                break;
            }
            case 2: //match same
            {
                best_match_len = buff_bm.decodeFast(&pack.huf_match_lens[ones_group]);
                
                if(decoded_bytes + best_match_len  <=  offset)
                {
                    decoded_bytes += best_match_len;
                }
                else if(decoded_bytes + best_match_len <= offset + length)
                {
                    buff_bm.getBuffer(bit0, left);
                    curr_pos = pack.bm.GetPos() - 1;
                    bits_read = 8 - pack.bm.GetWordPos();
                    nesting++;
                    decomp_vec_rrr_range(prev_vec_match, offset, decoded_bytes + best_match_len - offset, pos, decomp_data, start_id, true);
                    nesting--;
                    pack.bm.SetPos(curr_pos);
                    pack.bm.GetBitsAndDiscard(bits_read);
                    buff_bm.setBuffer(bit0, left);
                    
                    rest = decoded_bytes + best_match_len  - offset;
                    decoded_bytes = offset;
                    
                }
                else //if(decoded_bytes + best_match_len > offset + length)
                {
                    buff_bm.getBuffer(bit0, left);
                    
                    curr_pos = pack.bm.GetPos() - 1;
                    bits_read = 8 - pack.bm.GetWordPos();
                    nesting++;
                    decomp_vec_rrr_range(prev_vec_match, offset, length, pos, decomp_data, start_id, true);
                    nesting--;
                    pack.bm.SetPos(curr_pos);
                    pack.bm.GetBitsAndDiscard(bits_read);
                    buff_bm.setBuffer(bit0, left);
                    
                    decoded_bytes = offset;
                    rest = length;
                }
                
                break;
            }
            case 3: //zero run
            {
                zero_run_len = buff_bm.decodeFast(&(pack.huf_zeros_runs[ones_group]));
                if(zero_run_len > offset - decoded_bytes)
                {
                    rest = (zero_run_len - (offset - decoded_bytes)) ;
                    rest = rest < length ? rest : length;
                    zero_run_len = offset - decoded_bytes;
                }
                decoded_bytes += zero_run_len;
                if(rest)
                {
                    memcpy(decomp_data+pos, zeros_only_vector, rest);
                    pos = pos + rest;
                }
                break;
            }
            case 4: //ones run
            {
                
                ones_run_len = buff_bm.decodeFast(&(pack.huf_ones_runs[ones_group]));
                
                if(ones_run_len > offset - decoded_bytes)
                {
                    rest = (ones_run_len - (offset - decoded_bytes)) ;
                    rest = rest < length ? rest : length;
                    ones_run_len = offset - decoded_bytes;
                }
                
                decoded_bytes += ones_run_len;
                
                if(rest)
                {
                    memcpy(decomp_data+pos, ones_only_vector, rest);
                    pos = pos + rest;
                }
                break;
            }
            default: // Run of literals (2 (MIN_LITERAL_RUN) - MAX_LITERAL_RUN)
            {
                flag = flag - 3; // Shift by 3, to be able to difference between 2 and 3, 4 flags nad literal run of 2 or 3, 4
                
                if(flag < (int) (offset - decoded_bytes))
                {
                    
                    litRun = buff_bm.getBits(pack.used_bits_litRunSize[ones_group]);
                    if(!litRun) //litRun == 0 means used_bits_litRunSize bits were not enough to store size
                        litRun = buff_bm.getBits(pack.max_used_bits_litRunSize[ones_group]);
                    
                    litRun = litRun + pack.minLitRunSize[ones_group][flag] - 1; //1 is added to litRun, so it is never == 0
                    
                    // Skip run of literals
                    buff_bm.getBitsAndDiscard(litRun);
                    decoded_bytes += flag;
                }
                else
                {
                    // Discard description of size (bits) of run of literals
                    
                    litRun = buff_bm.getBits(pack.used_bits_litRunSize[ones_group]);
                    if(!litRun) //litRun == 0 means used_bits_litRunSize bits were not enough to store size
                        buff_bm.getBitsAndDiscard(pack.max_used_bits_litRunSize[ones_group]);
                    
                    rest = (flag - (offset - decoded_bytes)) ;
                    rest = rest < length ? rest : length;
                    
                    flag = offset - decoded_bytes;
                    for(int i = 0; i < flag; i++)
                    {
                        buff_bm.decodeFastLut(h_lit);
                    }
                    
                    for(int i = 0; i < (int) rest; i++)
                    {
                        decomp_data[pos++] =  buff_bm.decodeFastLut(h_lit);
                    }
                    decoded_bytes += flag + rest;
                }
                break;
            }
        }
        
    }
    decoded_bytes = rest;
    while(decoded_bytes < length)
    {
        flag = buff_bm.decodeFastLut(&pack.huf_flags);
        switch(flag)
        {
            case 0: //literal
            {
                decomp_data[pos++] =  buff_bm.decodeFastLut(h_lit);
                decoded_bytes++;
                break;
            }
            case 1: //match
            {
                // difference between current and match is stored
                tmp = buff_bm.decodeFastLut(&pack.huf_match_diff_MSB);
                best_pos = tmp << (pack.s.bit_size_id_match_pos_diff - MATCH_BITS_HUF);
                tmp = buff_bm.getBits(pack.s.bit_size_id_match_pos_diff - MATCH_BITS_HUF);
                best_pos = best_pos | tmp;
                
                best_pos += 1;
                
                best_pos = curr_non_copy_vec_id - best_pos;
            
                prev_vec_match = best_pos;
                best_match_len = buff_bm.decodeFast(&pack.huf_match_lens[ones_group]);
    
                if(best_match_len <= length - decoded_bytes)
                {
                    buff_bm.getBuffer(bit0, left);
                    curr_pos = pack.bm.GetPos() - 1;
                    bits_read = 8 - pack.bm.GetWordPos();
                    nesting++;
                    decomp_vec_rrr_range(best_pos, offset + decoded_bytes, best_match_len, pos, decomp_data, start_id, true);
                    nesting--;
                    pack.bm.SetPos(curr_pos);
                    pack.bm.GetBitsAndDiscard(bits_read);
                    decoded_bytes += best_match_len;
                    buff_bm.setBuffer(bit0, left);
                }
                else
                {
                    buff_bm.getBuffer(bit0, left);
                    
                    curr_pos = pack.bm.GetPos() - 1;
                    bits_read = 8 - pack.bm.GetWordPos();
                    nesting++;
                    decomp_vec_rrr_range(best_pos, offset + decoded_bytes, length - decoded_bytes, pos, decomp_data, start_id, true);
                    nesting--;
                    pack.bm.SetPos(curr_pos);
                    pack.bm.GetBitsAndDiscard(bits_read);
                    buff_bm.setBuffer(bit0, left);
                    decoded_bytes += length - decoded_bytes;
                }
                break;
            }
            case 2: //same match
            {
                best_match_len = buff_bm.decodeFast(&pack.huf_match_lens[ones_group]);
                
                if(best_match_len <= length - decoded_bytes)
                {
                    buff_bm.getBuffer(bit0, left);
                    
                    
                    curr_pos = pack.bm.GetPos() - 1;
                    bits_read = 8 - pack.bm.GetWordPos();
                    nesting++;
                    decomp_vec_rrr_range(prev_vec_match, offset + decoded_bytes, best_match_len, pos, decomp_data, start_id, true);
                    nesting--;
                    pack.bm.SetPos(curr_pos);
                    pack.bm.GetBitsAndDiscard(bits_read);
                    buff_bm.setBuffer(bit0, left);
                    decoded_bytes += best_match_len;
                }
                else
                {
                    buff_bm.getBuffer(bit0, left);
                    curr_pos = pack.bm.GetPos() - 1;
                    bits_read = 8 - pack.bm.GetWordPos();
                    nesting++;
                    decomp_vec_rrr_range(prev_vec_match, offset + decoded_bytes, length - decoded_bytes, pos, decomp_data, start_id, true);
                    nesting--;
                    pack.bm.SetPos(curr_pos);
                    pack.bm.GetBitsAndDiscard(bits_read);
                    buff_bm.setBuffer(bit0, left);
                    decoded_bytes += length - decoded_bytes;
                }
                break;
            }
            case 3: //zero run
            {
                
                zero_run_len = buff_bm.decodeFast(&(pack.huf_zeros_runs[ones_group]));
                zero_run_len = (zero_run_len <= length - decoded_bytes ) ? zero_run_len : length - decoded_bytes;
                memcpy(decomp_data + pos, zeros_only_vector, zero_run_len);
                decoded_bytes = decoded_bytes + zero_run_len;
                pos = pos + zero_run_len;
                break;
            }
            case 4: //one run
            {
                ones_run_len = buff_bm.decodeFast(&(pack.huf_ones_runs[ones_group]));
                ones_run_len = (ones_run_len <= length - decoded_bytes)  ? ones_run_len : length - decoded_bytes;
                memcpy(decomp_data + pos, ones_only_vector, ones_run_len);
                decoded_bytes = decoded_bytes + ones_run_len;
                pos = pos + ones_run_len;
                break;
                
            }
            default: // Run of literals (2-200)
            {
                flag = flag - 3; // Shift by 2, to be able to difference between 2 and 3 flags nad literal run of 2 or 3
                
                if(flag + decoded_bytes > length )
                    flag = length - decoded_bytes;
                
                litRun = buff_bm.getBits(pack.used_bits_litRunSize[ones_group]);
                if(!litRun) //litRun == 0 means used_bits_litRunSize bits were not enough to store size
                    buff_bm.getBitsAndDiscard(pack.max_used_bits_litRunSize[ones_group]);
                
                for(int i = 0; i < flag; i++)
                {
                    decomp_data[pos++] =  buff_bm.decodeFastLut(h_lit);
                }
                decoded_bytes += flag;
                break;
            }
        }
    }
    
    if(!is_unique_id && max_stored_unique)
    {
        uchar_t * vector;
        if(done_unique.size() > max_stored_unique)
        {
            toDelete = stored_unique.back();
            stored_unique.pop_back();
            vector = done_unique[toDelete];
            done_unique.erase(toDelete);
        }
        else
        {
            vector = new uchar_t[pack.s.vec_len];
        }
        
        memcpy(vector, decomp_data + (pos - pack.s.vec_len), pack.s.vec_len);
        done_unique[curr_non_copy_vec_id] = vector;
        stored_unique.push_front(curr_non_copy_vec_id);
    }
}

void Decompressor::decode_perm(int no_haplotypes, int vec2_start, uint32_t *perm, uchar_t *decomp_data_perm, uchar_t *decomp_data)
{
    for (uint32 j = 0; j < (uint32) no_haplotypes; ++j)
    {
        auto x = perm[j];
        uchar_t x_p = perm_lut[x % 8];
        uchar_t j_p = perm_lut[j % 8];
        int x8 = x / 8;
        int j8 = j / 8;
        
        if (decomp_data_perm[x8] & x_p)
            decomp_data[j8] += j_p;
        if (decomp_data_perm[vec2_start + x8] & x_p)
            decomp_data[vec2_start + j8] += j_p;
    }
}

void Decompressor::decode_perm_rev(int no_haplotypes, int vec2_start, uint32_t *rev_perm, uchar_t *decomp_data_perm, uchar_t *decomp_data)
{
#if 0
	for (uint32 x = 0; x < no_haplotypes; ++x)
    {
		auto j = rev_perm[x];
		uchar_t x_p = perm_lut[x % 8];
		uchar_t j_p = perm_lut[j % 8];
		int x8 = x / 8;
		int j8 = j / 8;
		if (decomp_data_perm[x8] & x_p)
			decomp_data[j8] += j_p;
		if (decomp_data_perm[vec2_start + x8] & x_p)
			decomp_data[vec2_start + j8] += j_p;
	}
	return;
#endif
	
    uint32 x;
    for (x = 0; x + 8 < (uint32) no_haplotypes;)
    {
        int x8 = x / 8;
        int d_x8 = decomp_data_perm[x8];
        int d2_x8 = decomp_data_perm[vec2_start + x8];

        if(!d_x8 && !d2_x8)
        {
            x += 8;
            continue;
        }
	
        int x_p = 1 << 7;
        for(int i = 0; i < 8; ++i, ++x, x_p >>= 1)
        {
            auto j = rev_perm[x];

            uchar_t j_p = perm_lut[j % 8];
            int j8 = j / 8;
            
            if (d_x8 & x_p)
                decomp_data[j8] += j_p;
            if (d2_x8 & x_p)
                decomp_data[vec2_start + j8] += j_p;
        }
    }
	
    int x8 = x / 8;
    int d_x8 = decomp_data_perm[x8];
    int d2_x8 = decomp_data_perm[vec2_start + x8];

    for (; x < (uint32) no_haplotypes; ++x)
    {
        auto j = rev_perm[x];
        uchar_t x_p = perm_lut[x % 8];
        uchar_t j_p = perm_lut[j % 8];

        int j8 = j / 8;
        if (d_x8 & x_p)
            decomp_data[j8] += j_p;
        if (d2_x8 & x_p)
            decomp_data[vec2_start + j8] += j_p;
    }
}

void inline Decompressor::reverse_perm(uint32_t *perm, uint32_t *rev_perm, int no_haplotypes)
{
    for(int i = 0; i < no_haplotypes; ++i)
               rev_perm[perm[i]] = i;
}

#ifdef DEVELOPMENT_MODE
int Decompressor::initBVOut(bool header_on)
{
    
    if(out_name != "")
    {
        if(out_type == BV)
        {
            char *gz_fname = (char*) malloc(strlen(out_name.c_str())+10);
            snprintf(gz_fname,strlen(out_name.c_str())+10,"%s.dcmp.bv", out_name.c_str());
            bv_out = fopen(gz_fname, "wb");
            
            if(header_on)
            {
                fwrite(&pack.s.n_samples, sizeof(uint32_t), 1, bv_out);
                fwrite(&pack.s.ploidy, sizeof(uint32_t), 1, bv_out);
                uint32_t no_var = pack.no_vec/2;
                fwrite(&no_var, sizeof(uint32_t), 1, bv_out);
            }
            free(gz_fname);
        }
        else if(out_type == TXT_BV)
        {
            char *gz_fname = (char*) malloc(strlen(out_name.c_str())+15);
            snprintf(gz_fname,strlen(out_name.c_str())+15,"%s.dcmp.bv.txt", out_name.c_str());
            bv_out = fopen(gz_fname, "w");
            free(gz_fname);
        }
    }
    else
    {
        bv_out = stdout;
    }
    
    return 0;
}

void Decompressor::getBV()
{
    // Permutation LUT
    for(int i = 0; i < 8; ++i)
        perm_lut[i] = 1 << (7 - i);
    
    uchar_t * decomp_data = nullptr;
    uchar_t * decomp_data_perm = nullptr;
    uint32_t * perm = nullptr;
    uint32_t * rev_perm = nullptr;
    uint32_t pos;
    
    uint32_t no_haplotypes = pack.s.n_samples * pack.s.ploidy;
    perm = new uint32_t[no_haplotypes];
    rev_perm = new uint32_t[no_haplotypes];
    
    done_unique.clear();
    uint32_t g, vec1_start, vec2_start;
    
    uint32_t  end;
    uint32_t i = 0;
    
    uint32_t block_id, prev_block_id = 0xFFFF;
    buff_bm.setBitMemory(&pack.bm);
    
    if(no_haplotypes & 7)//%8)
    {
        end = pack.s.vec_len - 1;
        g = no_haplotypes & 7;//%8;
    }
    else
    {
        end = pack.s.vec_len;
        g = 0;
    }
    
    decomp_data = new uchar_t[pack.s.vec_len*2];
    decomp_data_perm = new uchar_t[pack.s.vec_len*2];
    
    for(i = 0; i < pack.no_vec; )
    {
        pos = 0;
        // Permutations
        block_id = i/pack.s.max_no_vec_in_block; // Pair of vectors always in the same block (so it is enough to find first vector block
        if(block_id != prev_block_id)
        {
            // Get approproate permutations
            pack.getPermArray(block_id, perm);
            reverse_perm(perm, rev_perm, no_haplotypes);
            
        }
        
        fill_n(decomp_data, pack.s.vec_len*2, 0);
        vec1_start = 0;
        vec2_start = pack.s.vec_len;
        decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
        decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
        
        decode_perm_rev(no_haplotypes, vec2_start, rev_perm, decomp_data_perm, decomp_data);
        
        if(out_type == BV)
            fwrite(decomp_data, pack.s.vec_len*2, 1, bv_out);
        else //TXT_BV
        {
            for(uint32_t b = 0; b < pack.s.vec_len*2; b++)
                fprintf(bv_out, "%d ", (int)decomp_data[b]);
            fprintf(bv_out, "\n");
            
        }
        prev_block_id = block_id;
    }
    
    for (auto & it:done_unique)
    {
        delete [] it.second;
    }
    done_unique.clear();
    
    if(decomp_data)
        delete [] decomp_data;
    if(decomp_data_perm)
        delete [] decomp_data_perm;
    if(perm)
        delete [] perm;
    if(rev_perm)
        delete [] rev_perm;
}

void Decompressor::getVarBV(uint32_t var_to_dec)
{
    if(var_to_dec >= pack.no_vec/2)
    {
        cout << "Variant number "<< var_to_dec << " out of scope" << endl;
        exit(2);
    }

    // Permutation LUT
    for(int i = 0; i < 8; ++i)
        perm_lut[i] = 1 << (7 - i);
    
    uchar_t * decomp_data = nullptr;
    uchar_t * decomp_data_perm = nullptr;
    uint32_t * perm = nullptr;
    uint32_t * rev_perm = nullptr;
    uint32_t pos;
   
    
    uint32_t no_haplotypes = pack.s.n_samples * pack.s.ploidy;
    perm = new uint32_t[no_haplotypes];
    rev_perm = new uint32_t[no_haplotypes];
    
    done_unique.clear();
    uint32_t g, vec1_start, vec2_start;
    
    uint32_t  end;
    
    uint32_t block_id, prev_block_id = 0xFFFF;
    buff_bm.setBitMemory(&pack.bm);
    
    if(no_haplotypes & 7)//%8)
    {
        end = pack.s.vec_len - 1;
        g = no_haplotypes & 7;//%8;
    }
    else
    {
        end = pack.s.vec_len;
        g = 0;
    }
    
    decomp_data = new uchar_t[pack.s.vec_len*2];
    decomp_data_perm = new uchar_t[pack.s.vec_len*2];
    
    uint32_t i = var_to_dec* pack.s.ploidy;
    {
        pos = 0;
        // Permutations
        block_id = i/pack.s.max_no_vec_in_block; // Pair of vectors always in the same block (so it is enough to find first vector block
        if(block_id != prev_block_id)
        {
            // Get approproate permutations
            pack.getPermArray(block_id, perm);
            reverse_perm(perm, rev_perm, no_haplotypes);
        }
        
        fill_n(decomp_data, pack.s.vec_len*2, 0);
        
        vec1_start = 0;
        vec2_start = pack.s.vec_len;
        decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
        decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
        
        decode_perm_rev(no_haplotypes, vec2_start, rev_perm, decomp_data_perm, decomp_data);
        
        if(out_type == BV)
            fwrite(decomp_data, pack.s.vec_len*2, 1, bv_out);
        else //TXT_BV
        {
            for(uint32_t b = 0; b < pack.s.vec_len*2; b++)
                fprintf(bv_out, "%d ", (int)decomp_data[b]);
            fprintf(bv_out, "\n");
            
        }
        prev_block_id = block_id;
    }
    
    for (auto & it:done_unique)
    {
        delete [] it.second;
    }
    done_unique.clear();
    
    if(decomp_data)
        delete [] decomp_data;
    if(decomp_data_perm)
        delete [] decomp_data_perm;
    if(perm)
        delete [] perm;
    if(rev_perm)
        delete [] rev_perm;
}

void Decompressor::getVarSampleBV(uint32_t var_to_dec, uint32_t sample_to_dec)
{
    if(var_to_dec >= pack.no_vec/2)
    {
        cout << "Variant number "<< var_to_dec << " out of scope" << endl;
        exit(2);
    }
    if(sample_to_dec * pack.s.ploidy >= pack.s.vec_len * 8)
    {
        cout << "Sample id "<< sample_to_dec << " out of scope" << endl;
        exit(2);
    }

    initialLut();
    
    uchar_t * decomp_data = nullptr;
    uchar_t * decomp_data_perm = nullptr;
    uint32_t * perm = nullptr;
    uint32_t * rev_perm = nullptr;
    uint32_t pos;
    
    uint32_t no_haplotypes = pack.s.n_samples * pack.s.ploidy;
    perm = new uint32_t[no_haplotypes];
    rev_perm = new uint32_t[no_haplotypes];
    
    done_unique.clear();
    uint32_t g, vec1_start, vec2_start;
    
    uint32_t  end;
    
    uint32_t block_id, prev_block_id = 0xFFFF;
    buff_bm.setBitMemory(&pack.bm);
    
    if(no_haplotypes & 7)//%8)
    {
        end = pack.s.vec_len - 1;
        g = no_haplotypes & 7;//%8;
    }
    else
    {
        end = pack.s.vec_len;
        g = 0;
    }
    
    decomp_data = new uchar_t[pack.s.vec_len*2];
    decomp_data_perm = new uchar_t[pack.s.vec_len*2];
    
    uint32_t i = var_to_dec* pack.s.ploidy;

    pos = 0;
    // Permutations
    block_id = i/pack.s.max_no_vec_in_block; // Pair of vectors always in the same block (so it is enough to find the first vector block)
    if(block_id != prev_block_id)
    {
        // Get approproate permutations
        pack.getPermArray(block_id, perm);
        reverse_perm(perm, rev_perm, no_haplotypes);
    }
    
    fill_n(decomp_data, pack.s.vec_len*2, 0);
   
    vec1_start = 0;
    vec2_start = pack.s.vec_len;
    decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
    decomp_vec_rrr_range(i++, 0, pack.s.vec_len, pos, decomp_data_perm, 0, false);
    
    decode_perm_rev(no_haplotypes, vec2_start, rev_perm, decomp_data_perm, decomp_data);
    
    for(int p= 0; p < pack.s.ploidy; p++)
    {
        vec1_start = (sample_to_dec*pack.s.ploidy + p)/8; //2 vectors per variant
        uchar gt = *(lut[decomp_data[vec1_start]][decomp_data[vec1_start + pack.s.vec_len]] + (sample_to_dec*pack.s.ploidy + p)%8);
        
        if(gt == bcf_gt_missing)
        {
            gt = '.';
            if(out_type == BV)
                fwrite(&gt, sizeof(uchar), 1, bv_out);
            else //TXT_BV
                fprintf(bv_out, ".");
        }
        else
        {
            gt = bcf_gt_allele(gt);
            if(out_type == BV)
                fwrite(&gt, sizeof(uchar), 1, bv_out);
            else //TXT_BV
                fprintf(bv_out, "%d", (int)gt);
        }
    }

    if(out_type == TXT_BV)
          fprintf(bv_out, "\n");
    
    prev_block_id = block_id;

    for (auto & it:done_unique)
    {
        delete [] it.second;
    }
    done_unique.clear();
    
    if(decomp_data)
        delete [] decomp_data;
    if(decomp_data_perm)
        delete [] decomp_data_perm;
    if(perm)
        delete [] perm;
    if(rev_perm)
        delete [] rev_perm;
}

void Decompressor::getSampleBV(uint32_t sample_to_dec)
{
    if(sample_to_dec * pack.s.ploidy >= pack.s.vec_len * 8)
    {
        cout << "Sample id "<< sample_to_dec << " out of scope" << endl;
        exit(2);
    }
    
    initialLut();
    uchar_t * res = nullptr, *resUnique = nullptr;
    res = new uchar_t[pack.no_vec];
    resUnique = new uchar_t[pack.no_non_copy];
    buff_bm.setBitMemory(&pack.bm);
    bool is_unique = false;
    uint32_t unique_pos = 0;
    uint64_t  curr_zeros = 0, curr_copy = 0;
    
    uint64_t bv_size = pack.rrr_zeros_only_bit_vector[0].size();
    
    zeros_only_bit_vector[0] = sdsl::bit_vector(bv_size);
    zeros_only_bit_vector[1] = sdsl::bit_vector(bv_size);
    copy_bit_vector[0] = sdsl::bit_vector(bv_size);
    copy_bit_vector[1] = sdsl::bit_vector(bv_size);
    
    uint64_t v_pos;
    for(v_pos = 0; v_pos + 64 < bv_size; v_pos += 64)
    {
        zeros_only_bit_vector[0].set_int(v_pos, pack.rrr_zeros_only_bit_vector[0].get_int(v_pos, 64), 64);
        zeros_only_bit_vector[1].set_int(v_pos, pack.rrr_zeros_only_bit_vector[1].get_int(v_pos, 64), 64);
        copy_bit_vector[0].set_int(v_pos, pack.rrr_copy_bit_vector[0].get_int(v_pos, 64), 64);
        copy_bit_vector[1].set_int(v_pos, pack.rrr_copy_bit_vector[1].get_int(v_pos, 64), 64);
    }
    
    uint64_t tail_len = bv_size - v_pos;
    zeros_only_bit_vector[0].set_int(v_pos, pack.rrr_zeros_only_bit_vector[0].get_int(v_pos, tail_len), tail_len);
    zeros_only_bit_vector[1].set_int(v_pos, pack.rrr_zeros_only_bit_vector[1].get_int(v_pos, tail_len), tail_len);
    copy_bit_vector[0].set_int(v_pos, pack.rrr_copy_bit_vector[0].get_int(v_pos, tail_len), tail_len);
    copy_bit_vector[1].set_int(v_pos, pack.rrr_copy_bit_vector[1].get_int(v_pos, tail_len), tail_len);
    
    uint32_t * perm = nullptr;
    
    perm = new uint32_t[pack.s.n_samples * pack.s.ploidy];
    
    uint32_t block_id, prev_block_id = 0xFFFF;
    
    uint32_t size = (pack.no_vec * pack.s.ploidy)/8 + ((pack.no_vec * pack.s.ploidy)%8 > 0 ? 1 : 0);
    CBitMemory dec_data;
    dec_data.Create(size);
    char posInByteMask = 0;
    char shift = 0;
    
    uint32_t ind_id = 0, ind_id_orig;

    for (uint32_t p = 0; p < pack.s.ploidy; p++ )
    {
        ind_id_orig = sample_to_dec * pack.s.ploidy + p;
        
        
        unique_pos = 0;
        curr_zeros = 0, curr_copy = 0;
        for (uint32_t i = 0; i < pack.no_vec; i++ )
        {
            // Permutations
            block_id = i/pack.s.max_no_vec_in_block; //pair of vectors always in the same block (so it is enough to find first vector block)
            if(block_id != prev_block_id)
            {
                //get approproate permutations
                pack.getPermArray(block_id, perm);
                
                ind_id = perm[ind_id_orig];
                
                posInByteMask = 1<<(7 - (ind_id)%8);
                shift = 7 - (ind_id)%8;
            }
            
            res[i] = get_vec_byte(i, ind_id/8, resUnique, is_unique, curr_zeros, curr_copy);
            
            if(is_unique)
                resUnique[unique_pos++] = res[i];
            
            // Permutations
            prev_block_id = block_id;
            dec_data.PutBit((res[i] & posInByteMask) >> shift);
        }
    }
    dec_data.FlushPartialWordBuffer();

    if(out_type == BV)
        fwrite(dec_data.mem_buffer, dec_data.mem_buffer_pos, 1, bv_out);
    else //TXT_BV
    {
        for(uint32_t b = 0; b <   dec_data.mem_buffer_pos; b++)
        {
            fprintf(bv_out, "%d ", (int)dec_data.mem_buffer[b]);
            if(b == pack.no_vec/8)
                fprintf(bv_out, "\n");
        }
        fprintf(bv_out, "\n");
    }
   
    if(res)
        delete [] res;
    if(resUnique)
        delete [] resUnique;
    if(perm)
        delete [] perm;
   
}


/************************/
uchar_t Decompressor::get_vec_byte(uint32 vec_id, uint32 byte_no, uchar_t * resUnique, bool & is_uniqe_id, uint64_t & curr_zeros, uint64_t & curr_copy)
{
    uint32_t tmp = 0;
    
    uint8_t parity = vec_id&1;
    uint64_t vector = vec_id>>1, curr_non_copy_vec_id;
    
    if(parity)
    {
        if(copy_bit_vector[parity][vector])
        {
            unsigned long long bit_pos = (curr_copy)*pack.used_bits_cp;
            
            pack.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
            pack.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos&7); // %8
            pack.bm_comp_copy_orgl_id.GetBits(tmp, pack.used_bits_cp);
            
            // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
            curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy;
            
            curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;
            is_uniqe_id = false;
            curr_copy++;
            return resUnique[curr_non_copy_vec_id];
        }
        else if(zeros_only_bit_vector[1][vector])
        {
            is_uniqe_id = false;
            curr_zeros++;
            return 0;
        }
        else
        {
            curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy;
            is_uniqe_id = true;
        }
    }
    else if(!zeros_only_bit_vector[0][vector])
    {
        is_uniqe_id = false;
        curr_zeros++;
        return 0;
    }
    else
    {
        if(copy_bit_vector[parity][vector])
        {
            unsigned long long bit_pos = (curr_copy)*pack.used_bits_cp;
            
            pack.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
            pack.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos&7); // %8
            pack.bm_comp_copy_orgl_id.GetBits(tmp, pack.used_bits_cp);
            
            // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
            curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy;
            
            curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;
            is_uniqe_id = false;
            curr_copy++;
            return resUnique[curr_non_copy_vec_id];
        }
        else
        {
            curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy;
            is_uniqe_id = true;
        }
    }
    
    uint32_t curr_pos;
    tmp = 0;
    unsigned long long full_pos = curr_non_copy_vec_id/FULL_POS_STEP  * (sizeof(uint32_t) + pack.used_bits_noncp*(FULL_POS_STEP-1)/BITS_IN_BYTE);
    pack.bm_comp_pos.SetPos(full_pos);
    pack.bm_comp_pos.GetWord(curr_pos);
    uint32_t j = ((curr_non_copy_vec_id-1)%FULL_POS_STEP)/(sizeof(uint32_t)*BITS_IN_BYTE), end = curr_non_copy_vec_id%FULL_POS_STEP;
    if(j)
    {
        pack.bm_comp_pos.SetPos(full_pos + sizeof(uint32_t) + j*pack.used_bits_noncp*sizeof(uint32_t));
        j = j*(sizeof(uint32_t)*BITS_IN_BYTE);
        
        if(end > j + 1)
            pack.bm_comp_pos.GetBitsAndDiscard((int32_t)pack.used_bits_noncp * (end - 1 - j));
        if(end  > j)
            pack.bm_comp_pos.GetBits(tmp, (int32_t)pack.used_bits_noncp);
    }
    else
    {
        if(end > 1)
            pack.bm_comp_pos.GetBitsAndDiscard((int32_t)pack.used_bits_noncp*(end - 1));
        if(end)
            pack.bm_comp_pos.GetBits(tmp, (int32_t)pack.used_bits_noncp);
        
    }
    curr_pos += tmp;
    
    uint32_t ones_group;
    uint32_t decoded_bytes = 0;
    uint32_t zero_run_len,ones_run_len;
    uint32_t best_pos = 0;
    uint32_t best_match_len = 0;
    int32_t  flag;
    uint32 litRun;
    
    buff_bm.SetPos(curr_pos);
    
    ones_group = buff_bm.decodeFastLut(&pack.huf_group_type);
    
    CHuffman * h_lit = &pack.huf_literals[ones_group];
    
    while(decoded_bytes < byte_no)
    {
        flag = buff_bm.decodeFastLut(&pack.huf_flags);
        
        switch(flag)
        {
            case 0: // literal
            {
                buff_bm.decodeFastLut(h_lit);
                decoded_bytes++;
                break;
            }
            case 1: // match
            {
                // Difference between current and match is stored
                tmp = buff_bm.decodeFastLut(&pack.huf_match_diff_MSB);
                best_pos = tmp << (pack.s.bit_size_id_match_pos_diff - MATCH_BITS_HUF);
                tmp = buff_bm.getBits(pack.s.bit_size_id_match_pos_diff - MATCH_BITS_HUF);
                best_pos = best_pos | tmp;
                
                best_pos += 1;
                best_pos = curr_non_copy_vec_id - best_pos;
                
                best_match_len = buff_bm.decodeFast(&pack.huf_match_lens[ones_group]);
                
                if(decoded_bytes + best_match_len  <=  byte_no)
                {
                    decoded_bytes += best_match_len;
                }
                else
                    return  resUnique[best_pos];
                
                break;
            }
            case 2: // match same
            {
                best_match_len = buff_bm.decodeFast(&pack.huf_match_lens[ones_group]);
                if(decoded_bytes + best_match_len  <=  byte_no)
                {
                    decoded_bytes += best_match_len;
                }
                else
                    return  resUnique[best_pos];
                break;
            }
            case 3: // zero run
            {
                zero_run_len = buff_bm.decodeFast(&(pack.huf_zeros_runs[ones_group]));
                
                if(decoded_bytes + zero_run_len <= byte_no)
                {
                    decoded_bytes += zero_run_len;
                }
                else
                    return 0;
                
                
                break;
            }
            case 4: // one run
            {
                ones_run_len = buff_bm.decodeFast(&(pack.huf_ones_runs[ones_group]));
                
                if(decoded_bytes + ones_run_len <= byte_no)
                {
                    decoded_bytes += ones_run_len;
                }
                else
                    return 0xFF;
                
                
                break;
            }
            default: // Run of literals (2 (MIN_LITERAL_RUN) - MAX_LITERAL_RUN)
            {
                flag = flag - 3; //shift by 3, to be able to difference between 2,3 and 4 flags nad literal run of 2,3 or 4
                
                if(flag + decoded_bytes <= byte_no)
                {
                    litRun = buff_bm.getBits(pack.used_bits_litRunSize[ones_group]);
                    if(!litRun) //litRun == 0 means used_bits_litRunSize bits were not enough to store size
                        litRun = buff_bm.getBits(pack.max_used_bits_litRunSize[ones_group]);
                    
                    litRun = litRun + pack.minLitRunSize[ones_group][flag] - 1; //1 is added to litRun, so it is never == 0
                    
                    // Skip run of literals
                    buff_bm.getBitsAndDiscard(litRun);
                    decoded_bytes += flag;
                }
                else
                {
                    // Discard description of size (bits) of run of literals
                    litRun = buff_bm.getBits(pack.used_bits_litRunSize[ones_group]);
                    if(!litRun) //litRun == 0 means used_bits_litRunSize bits were not enough to store size
                        buff_bm.getBitsAndDiscard(pack.max_used_bits_litRunSize[ones_group]);
                    
                    flag = byte_no - decoded_bytes;
                    for(int i = 0; i < flag; i++)
                    {
                        buff_bm.decodeFastLut(h_lit);
                    }
                    return  buff_bm.decodeFastLut(h_lit);
                }
            }
        }
    }
    // Last byte
    {
        flag = buff_bm.decodeFastLut(&pack.huf_flags);
        
        switch(flag)
        {
            case 0:  //literal
            {
                return  buff_bm.decodeFastLut(h_lit);
                break;
            }
            case 1:  //match
            {
                tmp = buff_bm.decodeFastLut(&pack.huf_match_diff_MSB);
                best_pos = tmp << (pack.s.bit_size_id_match_pos_diff - MATCH_BITS_HUF);
                tmp = buff_bm.getBits(pack.s.bit_size_id_match_pos_diff - MATCH_BITS_HUF);
                best_pos = best_pos | tmp;
                
                best_pos += 1;
                best_pos = curr_non_copy_vec_id - best_pos;
                
                return  resUnique[best_pos];
                break;
            }
            case 2:  //same match
            {
                return  resUnique[best_pos];
                break;
            }
            case 3:
            {
                return 0;
                break;
            }
            case 4:
            {
                return 0xFF;
                break;
            }
            default:
            {
                // Discard description of size (bits) of run of literals
                litRun = buff_bm.getBits(pack.used_bits_litRunSize[ones_group]);
                if(!litRun) //litRun == 0 means used_bits_litRunSize bits were not enough to store size
                    buff_bm.getBitsAndDiscard(pack.max_used_bits_litRunSize[ones_group]);
                return  buff_bm.decodeFastLut(h_lit);
                break;
            }
        }
    }
    return 0;
}

#endif
