/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include "end_compressor.h"
#include <iostream>

void EndCompressor::Encode()
{
    
    rank_copy_bit_vector[0] = sdsl::rank_support_v5<>(&copy_bit_vector[0]);
    rank_copy_bit_vector[1] = sdsl::rank_support_v5<>(&copy_bit_vector[1]);
    
    rank_zeros_only_vector[0] = sdsl::rank_support_v5<>(&zeros_only_bit_vector[0]);
    rank_zeros_only_vector[1] = sdsl::rank_support_v5<>(&zeros_only_bit_vector[1]);
    
    // Calculate data stats (all histograms)
    calcModel();
   
    huf_literals = new CHuffman[s->ones_ranges];
    huf_zeros_runs = new CHuffman[s->ones_ranges];
    huf_ones_runs = new CHuffman[s->ones_ranges];
    huf_match_lens = new CHuffman[s->ones_ranges];
    
    uint32_t max_match_len = 1 << s->bit_size_match_len;
    
    for(int o_g = 0; o_g < (int) s->ones_ranges; o_g++)
    {
        // Huffman for literals
        huf_literals[o_g].Restart(256);
        for(uint64_t i = 0; i < 256; ++i)
        {
            if(hist_literals[o_g][i] <= UINT32_MAX)
                huf_literals[o_g].Insert((uint32)hist_literals[o_g][i]);
            else
            {
                cout << "Verify Huffman (variable to store frequency not large enough for literals)";
                cout << i << " " << hist_literals[o_g][i] << endl;
                exit(1);
            }
        }
        huf_literals[o_g].Complete();
        
        // Huffman for zeros runs
        uint32_t max_run_len = 1 << s->bit_size_run_len;
        huf_zeros_runs[o_g].Restart(max_run_len);
        for(uint64_t i = 0; i < max_run_len; ++i)
        {
            if((hist_zero_runs[o_g])[i] <= UINT32_MAX)
                huf_zeros_runs[o_g].Insert((uint32)hist_zero_runs[o_g][i]);
            else
            {
                cout << "Verify Huffman (variable to store frequency not large enough for zero runs)";
                cout << i << " " << hist_zero_runs[o_g][i] << endl;
                exit(1);
            }
        }
        huf_zeros_runs[o_g].Complete();
        
        // Huffman for ones runs
        uint32_t max_run_len_ones = 1 << s->bit_size_run_len;
        huf_ones_runs[o_g].Restart(max_run_len_ones);
        for(uint64_t i = 0; i < max_run_len_ones; ++i)
        {
            if((hist_ones_runs[o_g])[i] <= UINT32_MAX)
                huf_ones_runs[o_g].Insert((uint32)hist_ones_runs[o_g][i]);
            else
            {
                cout << "Verify Huffman (variable to store frequency not large enough for ones runs)";
                cout << i << " " << hist_ones_runs[o_g][i] << endl;
                exit(1);
            }
        }
        huf_ones_runs[o_g].Complete();
        
        // Huffman for match lengths
        huf_match_lens[o_g].Restart(max_match_len);
        for(uint64_t i = 0; i < max_match_len; ++i)
        {
            if(hist_match_lens[o_g][i] <= UINT32_MAX)
                huf_match_lens[o_g].Insert((uint32)hist_match_lens[o_g][i]);
            else
            {
                cout << "Verify Huffman (variable to store frequency not large enough for match lengths)";
                cout << i << " " << hist_match_lens[o_g][i] << endl;
                exit(1);
            }
        }
        huf_match_lens[o_g].Complete();
    }
    
    // Huffman for flags
    //shift by 3, to be able to difference between 2 and 3, 4 flags nad literal run of 2 or 3, 4
    uint64_t i;
    huf_flags.Restart(MAX_LITERAL_RUN+4);
    for(i = 0; i < MAX_LITERAL_RUN+4; ++i)
    {
        if((hist_flags[i]) <= UINT32_MAX)
            huf_flags.Insert((uint32)(hist_flags)[i]);
        else
        {
            cout << "Verify Huffman (variable to store frequency not large enough for flags)";
            cout << i << " " << hist_flags[i] << endl;
            exit(1);
        }
    }
    huf_flags.Complete();
    
    // Huffman for group_type
    huf_group_type.Restart(s->ones_ranges);
    for(i = 0; i < s->ones_ranges; ++i)
    {
        if((hist_group_type[i]) <= UINT32_MAX)
            huf_group_type.Insert((uint32)(hist_group_type)[i]);
        else
        {
            cout << "Verify Huffman (variable to store frequency not large enough for group_type)";
            cout << i << " " << hist_group_type[i] << endl;
            exit(1);
        }
    }
    huf_group_type.Complete();
    
    bm_huff.Create(bm.mem_buffer_pos);
    
    bm.Restart();
    
    litRunSize = new uint16_t[literalRunCount];
    literalRunCount = 0;
    
    for(int o = 0; o < (int) s->ones_ranges; o++)
        for(int i = 0; i <= (int) MAX_LITERAL_RUN + 3; i++) //shift by 3, to be able to difference between 3 and 4 flags nad literal run of 3 or 4
            
            minLitRunSize[o][i] = INT32_MAX;
    
    block = 0;
    uint64_t unique = 0;
    
    // Check and store literal flag lenghts; calculate minLitRunSize for each group and flag
    for(uint64_t cur_vec_id = 0; cur_vec_id < unique_no; ++cur_vec_id)
    {
        if(cur_vec_id > unique + unique_in_block[block])
        {
            unique += unique_in_block[block];
            block++;
        }
        getLitRunSizes(cur_vec_id);
    }
    
    // Huffman for group_type match_diff_MSB
    huf_match_diff_MSB.Restart(1 << MATCH_BITS_HUF);
    for(i = 0; i < (1 << MATCH_BITS_HUF); ++i)
    {
        if((hist_match_diff_MSB[i]) <= UINT32_MAX)
            huf_match_diff_MSB.Insert((uint32)(hist_match_diff_MSB)[i]);
        else
        {
            cout << "Verify Huffman (variable to store frequency not large enough for hist_match_diff_MSB)";
            cout << i << " " << hist_match_diff_MSB[i] << endl;
            exit(1);
        }
    }
    huf_match_diff_MSB.Complete();
    
    // Decrease by minLitRunSize for each group and flag; get litRunSizeMax, get hist_litRunSize_usedBits
    literalRunCount = 0;
    bm.Restart();
  
    decreaseLitRunSizesbyMin(unique_no);
    
    for(int o = 0; o < (int) s->ones_ranges; o++)
    {
        max_used_bits_litRunSize[o] = (int)bits_used(litRunSizeMax[o]);
    }
    
    uint32_t total, total_min = INT32_MAX;
    
    for(int o = 0; o < (int) s->ones_ranges; o++)
    {
        total_min = INT32_MAX;
        for(uint32_t b = max_used_bits_litRunSize[o]; b > 0; --b)
        {
            total = 0;
            for(int a = 0; a <= (int) b; a++)
                total += hist_litRunSize_usedBits[o][a]*b;
            for(int a = b+1; a <= (int) max_used_bits_litRunSize[o]; a++)
                total += (max_used_bits_litRunSize[o]+b) * hist_litRunSize_usedBits[o][a];
            
            if(total < total_min)
            {
                total_min = total;
                used_bits_litRunSize[o] = b;
            }
        }
    }
    block = 0;
    unique = 0;
    literalRunCount = 0;
    bm.Restart();
    for(uint64_t cur_vec_id = 0; cur_vec_id < unique_no; ++cur_vec_id)
    {
        if(cur_vec_id == unique + unique_in_block[block])
        {
            unique += unique_in_block[block];
            block++;
        }
        encode_vec(cur_vec_id);
    }
    
    bm.Close();
    
    bm_huff.PutBits(0, MAX_HUF_LUT_LEN); // Padding, so SpeedupLUT for Huffman can work nice at the end
    bm_huff.FlushPartialWordBuffer();
    bm_comp_copy_orgl_id.Create(copy_no*4);
    uint64_t end, j;
  
    for(i = 0; i < copy_no; i++)
    {
        bm_comp_copy_orgl_id.PutBits(comp_pos_copy[i], (int32_t)used_bits_cp);
    }
    bm_comp_copy_orgl_id.FlushPartialWordBuffer();
    
    delete [] comp_pos_copy;
    comp_pos_copy = NULL;
    
    bm_comp_pos.Create(unique_no*4);
    used_bits_noncp = bits_used(max_pos_diff);
    for(i = 0; i < unique_no; i += FULL_POS_STEP)
    {
        bm_comp_pos.PutWord(comp_pos_non_copy[i]);
        end = i + FULL_POS_STEP < unique_no ? i + FULL_POS_STEP : unique_no;
        for(j = i + 1; j < end; j++)
        {
            bm_comp_pos.PutBits(comp_pos_non_copy[j], (int32_t)used_bits_noncp);
        }
        bm_comp_pos.FlushPartialWordBuffer();
    }
    
    delete [] comp_pos_non_copy;
    comp_pos_non_copy = NULL;
    
    for(int i = 0; i < (int) zeros_only_bit_vector[0].size(); i++)
        zeros_only_bit_vector[0][i] = !zeros_only_bit_vector[0][i];
    
    sdsl::util::clear(rank_zeros_only_vector[0]);
    sdsl::util::clear(rank_zeros_only_vector[1]);
    sdsl::util::clear(rank_copy_bit_vector[0]);
    sdsl::util::clear(rank_copy_bit_vector[1]);
    sdsl::util::clear(rank_unique);
}

void EndCompressor::AddBlock(int &id_block, unsigned char *compressed_block, size_t n_recs, size_t compressed_size, std::vector<int> &perm, std::vector<bool> &zeros, std::vector<bool> &copies, uint32_t * origin_of_copy)
{
    assert((uint32) id_block == no_blocks);
    no_blocks++;
    perms.push_back(perm);
    bm.PutBytes(compressed_block, compressed_size);
    
    uint64_t block_start_vec_id = curr_vec_id;
    uint64_t cur_copy_no = 0, cur_unique_no = 0, zzeros = 0;
    
    for(uint64_t i = 0; i < n_recs; i++)
    {
        zeros_only_bit_vector[curr_vec_id%2][curr_vec_id/2] = zeros[i];
        copy_bit_vector[curr_vec_id%2][curr_vec_id/2] = copies[i];
        
        if(copies[i])
            origin_of_copy[cur_copy_no++] += block_start_vec_id;
        else if(!zeros[i])
            cur_unique_no++;
        else
            zzeros++;
        
        curr_vec_id++;
    }
    memcpy(comp_pos_copy + copy_no, origin_of_copy, cur_copy_no*sizeof(uint32_t));
    
    copy_no += cur_copy_no;
    unique_no += cur_unique_no;
    unique_in_block.push_back(cur_unique_no);
    
    curr_pos += compressed_size;
}

void EndCompressor::encode_vec(uint64_t vec_id)
{
    
    uint32_t decoded_bytes;
    uint32_t zero_run_len, ones_run_len;
    uint32_t best_pos = 0;
    uint32_t best_match_len = 0, pos_diff;
    uint64_t tmp;
    
    uint32_t  byte, flag;
    
    bm.FlushInputWordBuffer();
    
    if(vec_id%FULL_POS_STEP == 0)
    {
        comp_pos_non_copy[vec_id] = bm_huff.GetPos();
    }
    else
    {
        pos_diff = bm_huff.GetPos() - comp_pos_non_copy[vec_id/FULL_POS_STEP * FULL_POS_STEP];
        comp_pos_non_copy[vec_id] = pos_diff;
        if(pos_diff > max_pos_diff)
            max_pos_diff = pos_diff;
    }
    
    uint32_t ones_group;
    bm.GetBits(ones_group, s->bit_size_ones_goup);
    bm_huff.PutBits(huf_group_type.codes[ones_group].code, huf_group_type.codes[ones_group].len);
    
    decoded_bytes = 0;
    
    while(decoded_bytes < s->vec_len)
    {
        bm.GetBits(flag, 8);
        bm_huff.PutBits(huf_flags.codes[flag].code, huf_flags.codes[flag].len);
        
        switch(flag)
        {
            case 0: //literal x1
                bm.GetBits(byte, 8); //get byte
                bm_huff.PutBits(huf_literals[ones_group].codes[byte].code, huf_literals[ones_group].codes[byte].len);
                decoded_bytes++;
                break;
            case 1: //match
            {
                bm.GetBits(best_pos, s->bit_size_id);
                best_pos += block*s->max_no_vec_in_block;
                // Difference between unique id is of current vector and unique id of match vector
                best_pos = vec_id - (best_pos \
                                     - rank_copy_bit_vector[0](best_pos/2+(best_pos%2)) - rank_copy_bit_vector[1](best_pos/2)  \
                                     - rank_zeros_only_vector[0](best_pos/2+(best_pos%2)) - rank_zeros_only_vector[1](best_pos/2));
                
                best_pos -= 1; //shift to not waste 1 value
  
                tmp = (best_pos >> (s->bit_size_id_match_pos_diff - MATCH_BITS_HUF)) & bm.n_bit_mask[MATCH_BITS_HUF];
                bm_huff.PutBits(huf_match_diff_MSB.codes[tmp].code, huf_match_diff_MSB.codes[tmp].len);
                bm_huff.PutBits(best_pos & bm.n_bit_mask[s->bit_size_id_match_pos_diff - MATCH_BITS_HUF], (s->bit_size_id_match_pos_diff - MATCH_BITS_HUF));
                
                bm.GetBits(best_match_len, s->bit_size_match_len);
                
                // Huffman
                bm_huff.PutBits(huf_match_lens[ones_group].codes[best_match_len].code, huf_match_lens[ones_group].codes[best_match_len].len);
                
                decoded_bytes += best_match_len;
                break;
            }
            case 2: //match same
            {
                bm.GetBits(best_match_len, s->bit_size_match_len);
                // Huffman
                bm_huff.PutBits(huf_match_lens[ones_group].codes[best_match_len].code, huf_match_lens[ones_group].codes[best_match_len].len);
                decoded_bytes += best_match_len;
                break;
            }
            case 3: //zero run
            {
                bm.GetBits(zero_run_len, s->bit_size_run_len);
                // Huffman
                bm_huff.PutBits(huf_zeros_runs[ones_group].codes[zero_run_len].code, huf_zeros_runs[ones_group].codes[zero_run_len].len);
                
                for(int p = 0; p < (int) zero_run_len; p++)
                {
                    if(decoded_bytes < s->vec_len)
                    {
                        decoded_bytes++;
                    }
                    else
                        break;
                }
                break;
            }
            case 4: //one run
            {
                bm.GetBits(ones_run_len, s->bit_size_run_len);
                // Huffman
                bm_huff.PutBits(huf_ones_runs[ones_group].codes[ones_run_len].code, huf_ones_runs[ones_group].codes[ones_run_len].len);
                for(int p = 0; p < (int) ones_run_len; p++)
                {
                    if(decoded_bytes < s->vec_len)
                    {
                        decoded_bytes++;
                    }
                    else
                        break;
                }
                break;
            }
            default: //run of 2(MIN_LITERAL_RUN) - MAX_LITERAL_RUN literals
            {
                flag = flag - 3; //shift by 3, to be able to difference between 2, 3 and 4 flags nad literal run of 2, 3 or 4
                // Instead of always using "max_used_bits_litRunSize", use used_bits_litRunSize bits; if not enough (8 bits of 0x00 as a flag), use max_used_bits_litRunSize bits additionally (with full size)
                if(litRunSize[literalRunCount] <= bm.n_bit_mask[used_bits_litRunSize[ones_group]])
                    bm_huff.PutBits(litRunSize[literalRunCount++], used_bits_litRunSize[ones_group]);
                else
                {
                    bm_huff.PutBits(0, used_bits_litRunSize[ones_group]);
                    bm_huff.PutBits(litRunSize[literalRunCount++], max_used_bits_litRunSize[ones_group]);
                }
                for(int i = 0; i < (int) flag; i++)
                {
                    bm.GetBits(byte, 8); //get byte
                    bm_huff.PutBits(huf_literals[ones_group].codes[byte].code, huf_literals[ones_group].codes[byte].len);
                    decoded_bytes++;
                }
                break;
            }
        }
    }
    bm_huff.FlushPartialWordBuffer();
}

// Only to get sizes of literals run (in bits); rest is ommitted, bm is not altered, set difference for match pos
void EndCompressor::getLitRunSizes(uint64_t vec_id)
{
    uint32_t decoded_bytes;
    uint32_t zero_run_len, ones_run_len;
    uint32_t best_pos = 0;
    uint32_t best_match_len = 0;
    uint32_t  byte, flag;
    uint32_t ones_group;
    
    bm.FlushInputWordBuffer();
    bm.GetBits(ones_group, s->bit_size_ones_goup);
    decoded_bytes = 0;
    while(decoded_bytes < s->vec_len)
    {
        bm.GetBits(flag, 8);
        switch(flag)
        {
            case 0: //literal x1
                
                bm.GetBits(byte, 8); //get byte
                decoded_bytes++;
                break;
            case 1: //match
            {
                
                bm.GetBits(best_pos, s->bit_size_id);
                best_pos += block*s->max_no_vec_in_block; //threads!
    
                best_pos = vec_id - (best_pos \
                                     - rank_copy_bit_vector[0](best_pos/2+(best_pos%2)) - rank_copy_bit_vector[1](best_pos/2)  \
                                     - rank_zeros_only_vector[0](best_pos/2+(best_pos%2)) - rank_zeros_only_vector[1](best_pos/2));
                best_pos -= 1;
                
                hist_match_diff_MSB[(best_pos >> (s->bit_size_id_match_pos_diff - MATCH_BITS_HUF)) & bm.n_bit_mask[MATCH_BITS_HUF]]++;
                
                bm.GetBits(best_match_len, s->bit_size_match_len);
                decoded_bytes += best_match_len;
                break;
            }
            case 2: //match same
            {
                bm.GetBits(best_match_len, s->bit_size_match_len);
                decoded_bytes += best_match_len;
                break;
            }
            case 3: //zero run
            {
                bm.GetBits(zero_run_len, s->bit_size_run_len);
                for(int p = 0; p < (int) zero_run_len; p++)
                {
                    if(decoded_bytes < s->vec_len)
                    {
                        decoded_bytes++;
                    }
                    else
                        break;
                }
                break;
            }
            case 4: //one run
            {
                bm.GetBits(ones_run_len, s->bit_size_run_len);
                for(int p = 0; p < (int) ones_run_len; p++)
                {
                    if(decoded_bytes < s->vec_len)
                    {
                        decoded_bytes++;
                    }
                    else
                        break;
                }
                break;
            }
            default: //run of  2(MIN_LITERAL_RUN) - MAX_LITERAL_RUN literals
            {
                flag = flag - 3; //shift by 3, to be able to difference between 2, 3 and 4 flags nad literal run of 2, 3 or 4
                // Count literal run size (in bits)
                uint32_t lit_run_size = 0;
                for(int i = 0; i < (int) flag; i++)
                {
                    bm.GetBits(byte, 8); //get byte
                    decoded_bytes++;
                    lit_run_size += huf_literals[ones_group].codes[byte].len;
                    
                }
                // Remember literals description size (in bits)
                litRunSize[literalRunCount++] = lit_run_size;
                if(lit_run_size < minLitRunSize[ones_group][flag])
                    minLitRunSize[ones_group][flag] = lit_run_size;
                
                break;
            }
        }
    }
}

void EndCompressor::decreaseLitRunSizesbyMin(uint64_t unique_no)
{
    for(uint64_t cur_vec_id = 0; cur_vec_id < unique_no; ++cur_vec_id)
    {
        uint32_t decoded_bytes;
        uint32_t zero_run_len, ones_run_len;
        uint32_t best_pos = 0;
        uint32_t best_match_len = 0;
        uint32_t  byte, flag;
        uint32_t ones_group;
        
        bm.FlushInputWordBuffer();
        bm.GetBits(ones_group, s->bit_size_ones_goup);
        decoded_bytes = 0;
        while(decoded_bytes < s->vec_len)
        {
            bm.GetBits(flag, 8);
            switch(flag)
            {
                case 0: //literal x1
                {
                    bm.GetBits(byte, 8); //get byte
                    decoded_bytes++;
                    break;
                }
                case 1: //match
                {
                    bm.GetBits(best_pos, s->bit_size_id);
                    bm.GetBits(best_match_len, s->bit_size_match_len);
                    decoded_bytes += best_match_len;
                    break;
                }
                case 2: //match same
                {
                    bm.GetBits(best_match_len, s->bit_size_match_len);
                    decoded_bytes += best_match_len;
                    break;
                }
                case 3: //zero run
                {
                    bm.GetBits(zero_run_len, s->bit_size_run_len);
                    for(int p = 0; p < (int) zero_run_len; p++)
                    {
                        if(decoded_bytes < s->vec_len)
                        {
                            decoded_bytes++;
                        }
                        else
                            break;
                    }
                    break;
                }
                case 4: //one run
                {
                    bm.GetBits(ones_run_len, s->bit_size_run_len);
                    for(int p = 0; p < (int) ones_run_len; p++) 
                    {
                        if(decoded_bytes < s->vec_len)
                        {
                            decoded_bytes++;
                        }
                        else
                            break;
                    }
                    break;
                }
                default: // Run of 2(MIN_LITERAL_RUN) - MAX_LITERAL_RUN literals
                {
                    flag = flag - 3; // Shift by 3, to be able to difference between 2, 3 and 4 flags nad literal run of 2, 3 or 4
                    // Count literal run size (in bits)
                    for(int i = 0; i < (int) flag; i++)
                    {
                        bm.GetBits(byte, 8); //get byte
                        decoded_bytes++;
                    }
                    // Remember literals description size (in bits)
                    litRunSize[literalRunCount] = litRunSize[literalRunCount] -  minLitRunSize[ones_group][flag] + 1; //so it is never == 0
                
                    // Stats
                    if(litRunSize[literalRunCount] > litRunSizeMax[ones_group])
                        litRunSizeMax[ones_group]  = litRunSize[literalRunCount];
                    
                    hist_litRunSize_usedBits[ones_group][(int) bits_used(litRunSize[literalRunCount])]++;
                    literalRunCount++;
                    break;
                }
            }
        }
    }
}

char EndCompressor::bits_used(unsigned int n)
{
    char bits = 0;
    while(n)
    {
        n = n >> 1; bits++;
    }
    return bits;
}

int EndCompressor::storeArchive(const char * arch_name)
{
    char *fname = (char*) malloc(strlen(arch_name)+5);
    snprintf(fname, strlen(arch_name)+5,"%s.gtc", arch_name);
 
    sdsl::osfstream out(fname, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        if (sdsl::util::verbose) {
            std::cerr<<"ERROR: store_to_file not successful for: `" << fname << "`" << std::endl;
        }
        exit(1);
    }
    sdsl::rrr_vector<> rrr_bit_vector[5];

    for(int v = 0; v < 2; v++)
    {
        rrr_bit_vector[v] = sdsl::rrr_vector<>(zeros_only_bit_vector[v]);
        rrr_bit_vector[v].serialize(out);
        sdsl::util::clear(rrr_bit_vector[v]);
        
    }
    for(int v = 0; v < 2; v++)
    {
        rrr_bit_vector[v+2] = sdsl::rrr_vector<>(copy_bit_vector[v]);
        rrr_bit_vector[v+2].serialize(out);
        sdsl::util::clear(rrr_bit_vector[v+2]);
        
    }
    out.close();
    if (sdsl::util::verbose) {
        std::cerr<<"INFO: store_to_file: `"<<fname<<"`"<<std::endl;
    }

    FILE * comp  = fopen(fname, "ab");
    if (!comp) {
        
        std::cerr<<"ERROR: storing archive not successful for: `"<<fname<<"`"<<std::endl;
        exit(1);
    }
    free(fname);
    
    fwrite(&s->ones_ranges, sizeof(uchar), 1, comp);
    fwrite(&s->ploidy, sizeof(uchar), 1, comp);
    fwrite(&s->vec_len, sizeof(s->vec_len), 1, comp);
    
    uchar * mem = nullptr;
    uint32_t len_huf;
    
    // Huffman for huf_match_diff_MSB flag
    huf_match_diff_MSB.StoreTree(mem, len_huf);
    fwrite(&len_huf, sizeof(uint32_t), 1, comp);
    fwrite(mem, sizeof(uchar), len_huf, comp);
    if(mem)
    {
        delete [] mem;
        mem = nullptr;
    }
    
    // Huffman for group_type flag
    huf_group_type.StoreTree(mem, len_huf);
    fwrite(&len_huf, sizeof(uint32_t), 1, comp);
    fwrite(mem, sizeof(uchar), len_huf, comp);
    if(mem)
    {
        delete [] mem;
        mem = nullptr;
    }

    // Huffman for literals zero and one runs, and match lens in ONES_RANGES groups
    for (int o_g = 0; o_g < (int) s->ones_ranges; o_g++)
    {
        huf_literals[o_g].StoreTree(mem, len_huf);
        fwrite(&len_huf, sizeof(uint32_t), 1, comp);
        fwrite(mem, sizeof(uchar), len_huf, comp);
        if(mem)
        {
            delete [] mem;
            mem = nullptr;
        }

        huf_zeros_runs[o_g].StoreTree(mem, len_huf);
        fwrite(&len_huf, sizeof(uint32_t), 1, comp);
        fwrite(mem, sizeof(uchar), len_huf, comp);
        if(mem)
        {
            delete [] mem;
            mem = nullptr;
        }

        huf_ones_runs[o_g].StoreTree(mem, len_huf);
        fwrite(&len_huf, sizeof(uint32_t), 1, comp);
        fwrite(mem, sizeof(uchar), len_huf, comp);
        if(mem)
        {
            delete [] mem;
            mem = nullptr;
        }

        huf_match_lens[o_g].StoreTree(mem, len_huf);
        fwrite(&len_huf, sizeof(uint32_t), 1, comp);
        fwrite(mem, sizeof(uchar), len_huf, comp);
        if(mem)
        {
            delete [] mem;
            mem = nullptr;
        }
    }
    
    // Flags for runs of literals
    huf_flags.StoreTree(mem, len_huf);
    fwrite(&len_huf, sizeof(uint32_t), 1, comp);
    fwrite(mem, sizeof(uchar), len_huf, comp);
    if(mem)
    {
        delete [] mem;
        mem = nullptr;
    }

    fwrite(&no_vec, sizeof(no_vec), 1, comp);
    fwrite(&copy_no, sizeof(copy_no), 1, comp);
    fwrite(&used_bits_cp, sizeof(char), 1, comp);
    fwrite(&bm_comp_copy_orgl_id.mem_buffer_pos, sizeof(int32_t), 1, comp);
    fwrite(bm_comp_copy_orgl_id.mem_buffer, 1, bm_comp_copy_orgl_id.mem_buffer_pos, comp);
    fwrite(&unique_no, sizeof(unique_no), 1, comp);
    
    for(int o = 0; o < (int) s->ones_ranges; o++)
    {
        fwrite(&max_used_bits_litRunSize[o], sizeof(uchar_t), 1, comp);
        fwrite(&used_bits_litRunSize[o], sizeof(uchar_t), 1, comp);
        fwrite(minLitRunSize[o], sizeof(uint32_t), MAX_LITERAL_RUN + 4, comp);
    }
    
    fwrite(&used_bits_noncp, sizeof(char), 1, comp);
    fwrite(&s->bit_size_id_match_pos_diff, sizeof(char), 1, comp);
    fwrite(&bm_comp_pos.mem_buffer_pos, sizeof(bm_comp_pos.mem_buffer_pos), 1, comp);
    fwrite(bm_comp_pos.mem_buffer, 1, bm_comp_pos.mem_buffer_pos, comp);
    
    // Permutations
    CBitMemory bv_perm;
    uint32_t bitsize_perm = bits_used(s->n_samples*s->ploidy);
    bv_perm.Create(((bitsize_perm*s->n_samples*s->ploidy)/8 + 1)* no_blocks);
    
    fwrite(&no_blocks, sizeof(no_blocks), 1, comp);
    fwrite(&s->max_no_vec_in_block, sizeof(s->max_no_vec_in_block), 1, comp);
    fwrite(&s->n_samples, sizeof(s->n_samples), 1, comp);
    
    for(uint32_t b = 0; b < no_blocks; ++b)
    {    for(uint32_t i = 0; i < s->n_samples*s->ploidy; ++i)
        {
            bv_perm.PutBits(perms[b][i], bitsize_perm);
            
        }
        bv_perm.FlushPartialWordBuffer();
    }
    fwrite(&bv_perm.mem_buffer_pos, 1, sizeof(bv_perm.mem_buffer_pos), comp);
    fwrite(bv_perm.mem_buffer, 1, bv_perm.mem_buffer_pos, comp);
    
    // Core (vector witch gt data)
    fwrite(bm_huff.mem_buffer, 1, bm_huff.mem_buffer_pos, comp);
    
    fclose(comp);
    no_vec= 0;
    bm_huff.Close();
    bm.Close();
    bm_comp_pos.Close();
    bm_comp_copy_orgl_id.Close();
    for (int o_g = 0; o_g < (int) s->ones_ranges; o_g++)
    {
        huf_literals[o_g].Restart();
        huf_zeros_runs[o_g].Restart();
        huf_ones_runs[o_g].Restart();
        huf_match_lens[o_g].Restart();
    }
    return 0;
}

void EndCompressor::calcModel()
{
    uint32_t comp_no_matches = 0;
    uint32_t comp_no_literals = 0;
    uint32_t comp_zero_run = 0;
    uint32_t comp_ones_run = 0;
    uint32_t n_copies = 0;
    uint32_t same_vec_match = 0;
    uint32_t comp_no_literals_in_run = 0;
    
    uint64_t i = 0, zeros_no  = 0;
    for ( i = 0; i < no_vec; i++)
    {
        if(zeros_only_bit_vector[i % 2][i / 2])
        {
            zeros_no++;
            continue;
        }
        else if(!copy_bit_vector[i% 2][i / 2])
        {
            unique[i] = 1;
        }
        else
            n_copies++;
    }
    rank_unique = sdsl::rank_support_v5<>(&unique);
    
    uint64_t curr_non_copy_vec_id = 0;
    uint64_t  copy_id = 0, origin_unique_id;
    uint32_t max_diff_copy = 0;
    for (uint64_t i = 0; i < no_vec; i++)
    {
        if(zeros_only_bit_vector[i % 2][i / 2])
            continue;
        if(!copy_bit_vector[i% 2][i / 2])
        {
            curr_non_copy_vec_id++;
            continue;
        }
        // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
        origin_unique_id = rank_unique(comp_pos_copy[copy_id]);
        
        // Store difference -1, to not waste one value
        comp_pos_copy[copy_id] = curr_non_copy_vec_id - origin_unique_id - 1;
        if(comp_pos_copy[copy_id] > max_diff_copy)
            max_diff_copy = comp_pos_copy[copy_id];
        copy_id++;
    }
    used_bits_cp = bits_used(max_diff_copy);
    
    // Histograms
    hist_match_diff_MSB = new uint64_t[1 << MATCH_BITS_HUF]();
    hist_zero_runs = new uint64_t*[s->ones_ranges];
    hist_ones_runs = new uint64_t*[s->ones_ranges];
    hist_literals = new uint64_t*[s->ones_ranges];
    hist_match_lens = new uint64_t*[s->ones_ranges];
    for(int o = 0; o < (int) s->ones_ranges; o++)
    {
        hist_literals[o] = new uint64_t[256]();
        hist_zero_runs[o] = new uint64_t[1<<s->bit_size_run_len]();
        hist_ones_runs[o] = new uint64_t[1<<s->bit_size_run_len]();
        hist_match_lens[o] = new uint64_t[1<<s->bit_size_match_len]();
    }
    
    hist_group_type = new uint64_t[s->ones_ranges]();
    uint32_t decoded_bytes;
    uint32_t ones_group, flag, best_match_len, best_pos, literal, zero_run_len, ones_run_len;
    bm.Restart();
    for(uint64_t i = 0; i < unique_no; i++)
    {
        bm.FlushInputWordBuffer();
        bm.GetBits(ones_group, s->bit_size_ones_goup);
        hist_group_type[ones_group]++;
        decoded_bytes = 0;
        while(decoded_bytes < s->vec_len)
        {
            bm.GetBits(flag, 8);
            hist_flags[flag]++;
            switch(flag)
            {
                case 0: //literal
                {
                    comp_no_literals++;
                    bm.GetBits(literal, 8);
                    hist_literals[ones_group][literal]++;
                    decoded_bytes++;
                    break;
                }
                case 1: //match
                {
                    comp_no_matches++;
                    bm.GetBits(best_pos, s->bit_size_id);
                    bm.GetBits(best_match_len, s->bit_size_match_len);
                    hist_match_lens[ones_group][best_match_len]++;
                    decoded_bytes += best_match_len;
                    break;
                }
                case 2: //match same
                {
                    same_vec_match++;
                    bm.GetBits(best_match_len, s->bit_size_match_len);
                    hist_match_lens[ones_group][best_match_len]++;
                    decoded_bytes += best_match_len;
                    break;
                }
                case 3:  //zero run
                {
                    comp_zero_run++;
                    bm.GetBits(zero_run_len, s->bit_size_run_len);
                    hist_zero_runs[ones_group][zero_run_len]++;
                    decoded_bytes += zero_run_len;
                    break;
                }
                case 4:  //ones run
                {
                    comp_ones_run++;
                    bm.GetBits(ones_run_len, s->bit_size_run_len);
                    hist_ones_runs[ones_group][ones_run_len]++;
                    decoded_bytes += ones_run_len;
                    break;
                }
                default: //5 +, literal run
                {
                    comp_no_literals += flag-3;
                    comp_no_literals_in_run += flag-3;
                    literalRunCount++;
                    for(int l = 0; l < (int) flag-3; l++)
                    {
                        bm.GetBits(literal, 8);
                        hist_literals[ones_group][literal]++;
                    }
                    decoded_bytes += (flag-3);
                    break;
                }
            }
        }
    }
    
    cout << "\nStats:\nno_vec:\t" <<  no_vec  << "\nno_unique:\t " <<  unique_no  << "\nno_zero_vectors:\t" <<  zeros_no << "\nno_copies:\t " <<  n_copies << "\nno_zero_runs:\t" <<   comp_zero_run << "\nno_ones_runs:\t" <<   comp_ones_run<< "\nno_matches:\t"<<  comp_no_matches << "\nno_same_vec_match\t" <<    same_vec_match << "\nno_literals(all):\t" <<     comp_no_literals<<  "\nno_lit_runs:\t" << literalRunCount << "\nno_literals_in_runs:\t" << comp_no_literals_in_run <<   endl <<   endl;
}
