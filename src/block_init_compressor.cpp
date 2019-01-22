/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include "block_init_compressor.h"
#include <list>
#include <utility>

void BlockInitCompressor::SetBlock(uint64_t _cur_no_vec, uchar_t * _data)
{
    data = _data;
    cur_no_vec = _cur_no_vec;
}

void BlockInitCompressor::PermuteBlock(vector<int> & perm, bool permute)
{
    if(permute)
    {
        permute_range_vec(0, cur_no_vec, perm);
    }
    else
    {
        perm.clear();
        perm.resize(s->vec_len * 8, 0);
        for (int i = 0; i < (int) s->vec_len * 8; ++i)
            perm[i] = i;
    }
}


bool BlockInitCompressor::Compress(vector<bool> &zeros, vector<bool> &copies,  uchar * & compressedBlock, size_t & compressed_size, uint32_t *& origin_of_copy)
{
    zeros_only = &zeros;
    copy = &copies;

    if(!allocated)
        Allocate();
    else
        Clear();
    
    uint64_t i;
    
    bm.Create(s->vec_len*s->max_no_vec_in_block/10);
    comp_no_matches = 0;
    comp_no_literals = 0;
    comp_zero_run = 0;
    comp_ones_run = 0;
    n_copies = 0;
    same_vec_match = 0;
    
    // Generate the lookup table
    for (i = 0; i < 256; i++)
    {
        lookup_t_ones[i] = (i & 1) + lookup_t_ones[i / 2];
    }
    
    for( cur_vec_id = 0; cur_vec_id < cur_no_vec; ++cur_vec_id)
    {
        comp_vec(cur_vec_id);
    }
    bm.FlushPartialWordBuffer();
    bm.TakeOwnership();
    compressedBlock = bm.mem_buffer;
    compressed_size = bm.mem_buffer_pos;
    origin_of_copy = new uint32_t[no_copy];
    
    memcpy(origin_of_copy, comp_pos_copy, no_copy*sizeof(uint32_t));
    
    bm.Close();
    return true;
}

void BlockInitCompressor::Allocate()
{
    uint64_t i;
    
    comp_pos_non_copy = new uint32_t[s->max_no_vec_in_block]();
    comp_pos_copy = new uint32_t[s->max_no_vec_in_block]();
    
    counters = new uchar[s->max_no_vec_in_block*s->vec_len]();
    
    // Allocate memory for HT keys
    ht = new uint32_t*[s->vec_len];
    ht_size    = new uint32_t[s->vec_len];
    ht_fill    = new uint32_t[s->vec_len];
    ht_entries = new uint32_t[s->vec_len];
    
    // Allocate memory for HT - hash for whole vector (used while searching for copies)
    vec_hash_size = 2 * s->max_no_vec_in_block;
    while (vec_hash_size & (vec_hash_size - 1))		// round to nearest power of two
        vec_hash_size &= vec_hash_size - 1;
    vec_hash_size *= 2;
    
    vec_hash = new uint32_t[vec_hash_size];
    fill_n(vec_hash, vec_hash_size, EMPTY);
    
    for (i = 0; i < s->vec_len; ++i)
    {
        ht[i] = new uint32_t[HT_INIT_SIZE];
        fill_n(ht[i], HT_INIT_SIZE, EMPTY);
        ht_size[i] = HT_INIT_SIZE;
        ht_fill[i] = 0;
        ht_entries[i] = 0;
    }
    
    temp_vec = new uchar[s->vec_len];
    allocated = true;
}

void BlockInitCompressor::Clear()
{
    uint64_t i;

    std::fill(zeros_only->begin(), zeros_only->end(), false);
    std::fill(copy->begin(), copy->end(), false);
    
    std::fill(comp_pos_non_copy, comp_pos_non_copy+s->max_no_vec_in_block, 0);
    std::fill(comp_pos_copy, comp_pos_copy+s->max_no_vec_in_block, 0);
    std::fill(counters, counters+s->max_no_vec_in_block*s->vec_len, 0);
    
    std::fill_n(vec_hash, vec_hash_size, EMPTY);
    
    for (i = 0; i < s->vec_len; ++i)
    {
        if(ht_size[i] != HT_INIT_SIZE)
        {
            delete[] ht[i];
            ht[i] = new uint32_t[HT_INIT_SIZE];
            ht_size[i] = HT_INIT_SIZE;
        }

        fill_n(ht[i], HT_INIT_SIZE, EMPTY);
        ht_fill[i] = 0;
        ht_entries[i] = 0;
    }
    
    no_non_copy = 0;
    no_copy = 0;
    no_zeros_only = 0;
    uniq_vec_counter = 0;
    n_vec_in_ht_parts = 0;
    n_vec_in_ht_vecs  = 0;
    
    idx_oldest_vec_in_ht_vecs = 0;
    idx_oldest_vec_in_ht_parts = 0;
}

void BlockInitCompressor::Deallocate()
{
    uint64_t i;

    if(allocated)
    {
        delete [] comp_pos_non_copy;
        delete [] comp_pos_copy;
        for (i = 0; i < s->vec_len; ++i)
        {
            delete[] ht[i];
        }
        delete [] ht;
        delete [] ht_size;
        delete [] ht_fill;
        delete [] ht_entries;
        delete [] vec_hash;
        delete [] temp_vec;
        delete [] counters;
    }
    allocated = false;
}

// Compression of single bit vector
void BlockInitCompressor::comp_vec(uint64_t vec_id)
{
    uint64_t i;

    uint64_t  curr_non_copy_vec_id;
    uint32_t pos_diff, lit_run_count = 0;
    uint64_t prev_vec_match = 1u << 30;
    
    // Search for vector copy
    uint64_t x = hash_fun(vec_id, 0, s->vec_len, vec_hash_size);
   
    if(x == ZEROS)
    {
        (*zeros_only)[vec_id] = 1;
        no_zeros_only++;
        return;
    }
    
    if(x != NO_POS)
        for (; vec_hash[x] != EMPTY; x = (x + 1) & (vec_hash_size - 1))
        {
            if (vec_hash[x] < idx_oldest_vec_in_ht_vecs)
                continue;
            if ((*zeros_only)[vec_hash[x]] || (*copy)[vec_hash[x]])
                continue;
            
            if (memcmp(&data[vec_id * s->vec_len], &data[vec_hash[x] * s->vec_len], s->vec_len) == 0)
            {
                // Vector is a copy of one of previous vector - nothing to encode, pointer to original vactor
                (*copy)[vec_id] = 1;
                comp_pos_copy[no_copy] = vec_hash[x];
                no_copy++;
                n_copies++;
                return;
            }
        }
    
    bm.FlushPartialWordBuffer();
    
    curr_non_copy_vec_id = vec_id - no_zeros_only - no_copy;
    if(curr_non_copy_vec_id%FULL_POS_STEP == 0)
    {
        comp_pos_non_copy[curr_non_copy_vec_id] = bm.GetPos();
    }
    else
    {
        pos_diff = bm.GetPos() - comp_pos_non_copy[curr_non_copy_vec_id/FULL_POS_STEP * FULL_POS_STEP];
        comp_pos_non_copy[curr_non_copy_vec_id] = pos_diff;
    }
    
    uchar_t ones_group = get_ones_group(vec_id);
    bm.PutBits(ones_group, s->bit_size_ones_goup);
    no_non_copy++;
    uniq_vec_counter++;

    if(x != NO_POS)
        ht_insert_vec(vec_id, x);
    
    fill_n(temp_vec, s->vec_len, 1);
    
    // Checking, if there is a zero run at the end - if so, it is encoded at the end
    uchar_t *data_all = &data[vec_id * s->vec_len];
    uint64_t zero_run_len_end = 0;
    for (i = s->vec_len - 1; zero_run_len_end < s->vec_len; i--, zero_run_len_end++)
        if (data_all[i])
            break;
    if(zero_run_len_end < MIN_ZERO_RUN_LEN)
        zero_run_len_end = 0;
    
    // Checking, if there is a one run at the end - if so, it is encoded at the end
    uint64_t ones_run_len_end = 0;
    for (i = s->vec_len - 1; ones_run_len_end < s->vec_len; i--, ones_run_len_end++)
        if (data_all[i] != 0xFF)
            break;
    if(ones_run_len_end < MIN_ONES_RUN_LEN)
        ones_run_len_end = 0;
    
    int64_t end1 = (HASH_KEY_LEN2 > zero_run_len_end ? HASH_KEY_LEN2 : zero_run_len_end);
    int64_t end2 = (HASH_KEY_LEN2 > ones_run_len_end ? HASH_KEY_LEN2 : ones_run_len_end);
    int64_t end = end1 > end2 ? end1 : end2;
    
    for(i = 0; (int64) (s->vec_len - i) >  end;)
    {
        uint64_t x;
        uint64_t best_pos = 0;
        uint64_t best_match_len = 0;
        uint64_t best_depth = 0;
        
        uchar_t *cur_data = &data[vec_id * s->vec_len + i];
        
        // Checking, if there is a zero run - if so, it is encoded
        uint32_t zero_run_len = 0;
        for (; i + zero_run_len < s->vec_len; zero_run_len++)
            if (cur_data[zero_run_len])
                break;

        if (zero_run_len >= MIN_ZERO_RUN_LEN)
        {
            if(lit_run_count >= MIN_LITERAL_RUN)
                encode_literal_run(lit_run_count, &data[vec_id * s->vec_len + i - lit_run_count]);
            lit_run_count = 0;
            
            bm.PutBits(3, 8);
            bm.PutBits(zero_run_len, s->bit_size_run_len);
            
            comp_zero_run++;
            fill_n(temp_vec + i, zero_run_len - MIN_ZERO_RUN_LEN, 0); // Indication, which positions should not be inserted in HT
            i += zero_run_len;
        
            continue;
        }
        
        // Checking, if there is a one run - if so, it is encoded
        uint32_t ones_run_len = 0;
        for (; i + ones_run_len < s->vec_len; ones_run_len++)
            if (cur_data[ones_run_len] != 0xFF)
                break;

        if (ones_run_len >= MIN_ONES_RUN_LEN)
        {
            if(lit_run_count >= MIN_LITERAL_RUN)
                encode_literal_run(lit_run_count, &data[vec_id * s->vec_len + i - lit_run_count]);
            lit_run_count = 0;
            
            bm.PutBits(4, 8);
            bm.PutBits(ones_run_len, s->bit_size_run_len);
    
            comp_ones_run++;
            fill_n(temp_vec + i, ones_run_len - MIN_ONES_RUN_LEN, 0);	// Indication, which positions should not be inserted in HT
            i += ones_run_len;
            continue;
        }
        
        if(s->vec_len - end <= 5) // Encoding iteral, if there is no point in looking for match
        {
            bm.PutBits(0, 8);
            bm.PutBits(data[vec_id * s->vec_len + i], 8);
            
            comp_no_literals++;
            counters[vec_id * s->vec_len + i] = 0;
            i++;
            lit_run_count++;
            
            if(lit_run_count == MAX_LITERAL_RUN)
            {
                encode_literal_run(lit_run_count, &data[vec_id * s->vec_len + i  - lit_run_count]);
                lit_run_count = 0;
            }
            continue;
        }
        
        // Searching for best match with HT
        if (!best_match_len)
        {
            x = hash_fun(vec_id, i, HASH_KEY_LEN2, ht_size[i]);
            if(x != NO_POS)
                for (; ht[i][x] != EMPTY; x = (x + 1) & (ht_size[i] - 1))
                {
                    if (ht[i][x] == REMOVED)
                        continue;
                    
                    uint64_t len;
                    uint32_t depth = 0;
                    uchar_t *vec_ptr = &data[ht[i][x] * s->vec_len + i];
                    uchar_t *counters_ptr = &counters[ht[i][x] * s->vec_len + i];
                    for (len = 0; len < s->vec_len - i; ++len)
                        if (cur_data[len] != vec_ptr[len])
                            break;
                        else if ((uint32_t)counters_ptr[len] > depth)
                        {
                            if ((uint32_t)counters_ptr[len] < s->max_depth)
                                depth = (uint32_t)counters_ptr[len];
                            else
                                break;
                        }
                    ht_checks++;
                    
                    if (len < HASH_KEY_LEN2)
                        continue;
                    
                    if ((len > best_match_len + 1) ||
                        (ht[i][x] == prev_vec_match && len + 1 >= best_match_len) ||
                        (best_pos != prev_vec_match && len > best_match_len) ||
                        (best_pos != prev_vec_match && len == best_match_len && depth < best_depth) ||
                        (best_pos != prev_vec_match && len == best_match_len && depth == best_depth && ht[i][x] > best_pos))
                    {
                        best_match_len = len;
                        best_depth = depth;
                        best_pos = ht[i][x];
                    }
                    
                    if (best_match_len + i >= s->vec_len)
                        break;
                }
        }
        
        if(best_match_len) //Encoding match
        {
            if (prev_vec_match == best_pos) //Encoding match to the same vector as previous match
            {
                same_vec_match++;
                if(lit_run_count >= MIN_LITERAL_RUN)
                    encode_literal_run(lit_run_count, &data[vec_id * s->vec_len + i - lit_run_count]);
                lit_run_count = 0;
            
                // match flag same
                bm.PutBits(2, 8);
                bm.PutBits((uint32_t)best_match_len, s->bit_size_match_len);
                
                for(uint64_t j = 0; j < best_match_len; ++j)
                    counters[vec_id * s->vec_len + i + j] = (uchar_t) (best_depth + 1);
                
                i += best_match_len;
            }
            else //Encoding match
            {
               
                if(lit_run_count >= MIN_LITERAL_RUN)
                    encode_literal_run(lit_run_count, &data[vec_id * s->vec_len + i - lit_run_count]);
                lit_run_count = 0;
                
                bm.PutBits(1, 8);
                bm.PutBits(best_pos, s->bit_size_id);
                bm.PutBits(best_match_len, s->bit_size_match_len);
                
                for(uint64_t j = 0; j < best_match_len; ++j)
                    counters[vec_id * s->vec_len + i + j] = (uchar_t) (best_depth + 1);
                
                comp_no_matches++;
                i += best_match_len;
            }
            prev_vec_match = best_pos;
        }
        else //Encoding literal
        {
            bm.PutBits(0, 8);
            lit_run_count++;
            bm.PutBits(data[vec_id * s->vec_len + i], 8);
            comp_no_literals++;
            counters[vec_id * s->vec_len + i] = 0;
            i++;
            
            if(lit_run_count == MAX_LITERAL_RUN)
            {
                encode_literal_run(lit_run_count, &data[vec_id * s->vec_len + i - lit_run_count]);
                lit_run_count = 0;
            }
        }
    }
    
    zero_run_len_end = zero_run_len_end > (s->vec_len - i)? s->vec_len - i : zero_run_len_end;
    ones_run_len_end = ones_run_len_end > (s->vec_len - i)? s->vec_len - i : ones_run_len_end;
    end = ones_run_len_end > zero_run_len_end ? ones_run_len_end : zero_run_len_end;

    // Encoding literals at the end of vector
    for(; i < s->vec_len - end; ++i)
    {
        bm.PutBits(0, 8); //literal
        bm.PutBits(data[vec_id * s->vec_len + i], 8);
        comp_no_literals++;
        counters[vec_id * s->vec_len + i] = 0;
        lit_run_count++;
        if(lit_run_count == MAX_LITERAL_RUN)
        {
            encode_literal_run(lit_run_count, &data[vec_id * s->vec_len + i + 1 - lit_run_count]);
            lit_run_count = 0;
        }
    }

    // Encoding literal island at the end
    if(lit_run_count >= MIN_LITERAL_RUN)
        encode_literal_run(lit_run_count, &data[vec_id * s->vec_len + i - lit_run_count]);
    lit_run_count = 0;
    
    if (zero_run_len_end && i <s->vec_len)	// Encoding zero run at the end of vector
    {
        bm.PutBits(3, 8);
        bm.PutBits(zero_run_len_end, s->bit_size_run_len);
        comp_zero_run++;
        fill_n(temp_vec + i, zero_run_len_end, 0);		// Indication, which positions should not be inserted in HT
    }
    
    if (ones_run_len_end && i <s->vec_len)	// Encoding one run at the end of vector
    {
        bm.PutBits(4, 8);
        bm.PutBits(ones_run_len_end, s->bit_size_run_len);
        comp_ones_run++;
        fill_n(temp_vec + i, ones_run_len_end, 0);		// Indication, which positions should not be inserted in HT
    }
    
    // Removing vector data from HT (only, if it was not a copy of other vector [only unique vectors])
    if (n_vec_in_ht_parts >= s->n_vec_history_parts)
    {
        while (true)
        {
            if (!(*zeros_only)[idx_oldest_vec_in_ht_parts] &&
                !(*copy)[idx_oldest_vec_in_ht_parts])
                break;
            ++idx_oldest_vec_in_ht_parts;
        }
        for (i = 0;  s->vec_len - i > HASH_KEY_LEN2 - 1; ++i)
            ht_remove(idx_oldest_vec_in_ht_parts, i);
        ++idx_oldest_vec_in_ht_parts;
    }
    else
        ++n_vec_in_ht_parts;
    
    // Store idx of oldest vector in HT
    if (n_vec_in_ht_vecs >= s->n_vec_history_vecs)
    {
        while (true)
        {
            if (!(*zeros_only)[idx_oldest_vec_in_ht_vecs] &&
                !(*copy)[idx_oldest_vec_in_ht_vecs])
                break;
            ++idx_oldest_vec_in_ht_vecs;
        }
        ++idx_oldest_vec_in_ht_vecs;
    }
    else
        ++n_vec_in_ht_vecs;
    
    // Inserting vector to HT (only, if it was not a copy of other vector [only unique vectors])
    for (i = 0;  s->vec_len - i > HASH_KEY_LEN2 - 1; ++i)
        if (temp_vec[i])
            ht_insert(vec_id, i);
}

void BlockInitCompressor::encode_literal_run(int32_t len, const uchar * lit_run)
{
    bm.discardBits(len*2*8);
    bm.PutBits(len + 3, 8);  //shift by 3, to be able to difference between 2, 3 and 4 flags nad literal run of 2, 3 or 4

    for(int i = 0; i < len; i++)
        bm.PutBits(lit_run[i], 8);
}

uchar_t BlockInitCompressor::get_ones_group(uint64_t vec_id)
{
    double count = 0;
    
    for (uint64_t i = 0; i < s->vec_len; ++i)
    {
        count += lookup_t_ones[data[vec_id * s->vec_len + i]];
    }
    
    if(count)
        count--;
    return (uchar_t)((count)/aux_dividor);
}

/************HT functions*************/

// Hashing functions; it calculates hash for position pointed by in_vec_pos
uint64_t BlockInitCompressor::hash_fun(uint64_t vec_id, uint64_t in_vec_pos, uint64_t len, uint64_t size)
{
    uint64_t r = 0;
 
    uchar_t *ptr_counters = &counters[vec_id * s->vec_len + in_vec_pos];
    uchar_t *ptr_data = &data[vec_id * s->vec_len + in_vec_pos];
    
    for (uint64_t i = 0; i < len; ++i)
    {
        if(ptr_counters[i] >= s->max_depth)
            return NO_POS;
        r = r * 0xcfcf + (uint64_t) ptr_data[i];
    }

    if(!r)
    {
        for (uint64_t i = 0; i < len; ++i)
        {
            if(ptr_data[i])
                return NO_POS;
        }
        return ZEROS;
    }

    return r & (size - 1);
}

bool BlockInitCompressor::ht_insert(uint64_t vec_id, uint64_t in_vec_pos)
{
    if (ht_fill[in_vec_pos] > ht_size[in_vec_pos] * 0.5)
        restruct(in_vec_pos);
    
    uint64_t idx = hash_fun(vec_id, in_vec_pos, HASH_KEY_LEN2, ht_size[in_vec_pos]);
    
    if (idx == NO_POS)
        return false;
    
    if (idx == ZEROS)
        return false;
    
    while (true)
    {
        if (ht[in_vec_pos][idx] == EMPTY || ht[in_vec_pos][idx] == REMOVED)
            break;
        idx++;
        if (idx == ht_size[in_vec_pos])
            idx = 0;
    }
    
    if (ht[in_vec_pos][idx] == EMPTY)
        ++ht_fill[in_vec_pos];
    ++ht_entries[in_vec_pos];
    
    ht[in_vec_pos][idx] = vec_id;
    
    return true;
}

// Removes old vector data from HT
bool BlockInitCompressor::ht_remove(uint64_t vec_id, uint64_t in_vec_pos)
{
    uint64_t idx = hash_fun(vec_id, in_vec_pos, HASH_KEY_LEN2, ht_size[in_vec_pos]);
    
    if (idx == NO_POS)
        return false;
    
    if (idx == ZEROS)
        return false;
    
    while (true)
    {
        if (ht[in_vec_pos][idx] == EMPTY || ht[in_vec_pos][idx] == vec_id)
            break;
        idx++;
        if (idx == ht_size[in_vec_pos])
            idx = 0;
    }
    
    if (ht[in_vec_pos][idx] == vec_id)
    {
        --ht_entries[in_vec_pos];
        ht[in_vec_pos][idx] = REMOVED;
    }
    
    return true;
}

// Inserting to hash table for whole vectors
// Resizing is not needed  - initial size is more than enough
bool BlockInitCompressor::ht_insert_vec(uint64_t vec_id, uint32_t key)
{
    while (vec_hash[key] != EMPTY)
        key = (key + 1) & (vec_hash_size - 1);
    
    vec_hash[key] = vec_id;
    return true;
}

void BlockInitCompressor::restruct(uint64_t in_vec_pos)
{
    delete[] ht[in_vec_pos];
    
    if (ht_entries[in_vec_pos] > ht_size[in_vec_pos] * 0.25)
        ht_size[in_vec_pos] *= 2;
    
    ht_fill[in_vec_pos] = 0;
    ht_entries[in_vec_pos] = 0;
    
    ht[in_vec_pos] = new uint32_t[ht_size[in_vec_pos]];
    fill_n(ht[in_vec_pos], ht_size[in_vec_pos], EMPTY);
    
    uint64_t start_vec_id = 0;
    if (cur_vec_id > s->n_vec_history_parts)
        start_vec_id = cur_vec_id - s->n_vec_history_parts;
    
    for (uint64_t i = start_vec_id; i < cur_vec_id; ++i)
        if((*zeros_only)[i] == 0 && (*copy)[i] == 0)
            ht_insert(i, in_vec_pos);
}

void BlockInitCompressor::permute_range_vec(uint64_t id_start, uint64_t id_stop, vector<int> &v_perm)
{
    int n_h_samples = v_perm.size();
    uint64_t part_vec = id_stop - id_start;
    vector<mc_vec_t> mc_vectors;
    mt19937 mt;
    
    // Calculate bit vectors with bits from the randomly chosen variants
    mc_vec_t empty_vec;
    empty_vec.fill(0);
    
    mc_vectors.resize(n_h_samples, empty_vec);
    
    vector<int> n_ones(n_h_samples, 0);
    vector<int> mc_ids;
	
    for(int i = 0; i < (int) part_vec; ++i)
    {
        uint32_t id_cur = i + id_start;
        auto cur_vec = data + id_cur * s->vec_len;
        bool empty = true;
            
        for (int j = 0; j < (int) s->vec_len; ++j)
            if (cur_vec[j])
            {
                empty = false;
                break;
            }
		
        if(!empty)
            mc_ids.push_back(i);		
    }
	
    random_shuffle(mc_ids.begin(), mc_ids.end());
    uint32_t mc_sort_size = mc_ids.size() > PART_TRIALS ? PART_TRIALS : mc_ids.size();
    sort(mc_ids.begin(), mc_ids.begin() + mc_sort_size);
	
    for (uint64_t i = 0; i < mc_sort_size; ++i)
    {
        uint32_t id_cur = mc_ids[i] + id_start;       
        
        auto cur_vec = data + id_cur * s->vec_len;
        int arr_id = i / 64;
        int arr_pos = i % 64;
        
        for (int j = 0; j < n_h_samples; ++j)
        {
            if(cur_vec[j / 8] & perm_lut8[j % 8])
            {
                mc_vectors[j][arr_id] += perm_lut64[arr_pos];
                ++n_ones[j];
            }
        }
    }
    
    // Calculate cost of the original permutation
    uint64_t cost_orig = 0;
    for (int i = 1; i < n_h_samples; ++i)
        cost_orig += bit_cost(mc_vectors[i - 1], mc_vectors[i]);
    
    // Determine the best permutation
    vector<int> perm;
    uint64_t cost_perm = 0;
    
    perm.push_back(0);
	
    list<pair<int, int>> density_list;

    density_list.push_back(make_pair(0, n_ones[0]));
    auto p = density_list.begin();
    for (int i = 1; i < (int) n_ones.size(); ++i)
        density_list.push_back(make_pair(i, n_ones[i]));

    // Insert two guards into list
    int huge_val = 1 << 28;
    density_list.push_back(make_pair(-1, -2*huge_val));
    density_list.push_back(make_pair(-1, 2*huge_val));

    density_list.sort([](pair<int, int> &x, pair<int, int> &y) {return x.second < y.second; });

    int n_cur_samples = n_h_samples;

    while (n_cur_samples)
    {
        uint64_t best_cost = huge_val;
        auto best_p = p;

        // Look for the most similar vector to *p
        // Starts from p and moves up and down on the list ordered according to the number of ones
        // This limits the number of vector pairs that must be evaluated
        auto p_down = p;
        auto p_up = p;
        --p_down;
        ++p_up;

        uint64_t dif_up = abs(p->second - p_up->second);
        uint64_t dif_down = abs(p->second - p_down->second);

        while (true)
        {
            uint64_t min_dif = min(dif_up, dif_down);

            if (min_dif >= best_cost)
                break;

            if (dif_up < dif_down)
            {
                uint64_t cost = bit_cost(mc_vectors[p->first], mc_vectors[p_up->first], best_cost);
                if (cost < best_cost)
                {
                    best_cost = cost;
                    best_p = p_up;

                    if (best_cost == 0)
			break;
                }

                ++p_up;
                dif_up = abs(p->second - p_up->second);
            }
            else
            {
                uint64_t cost = bit_cost(mc_vectors[p->first], mc_vectors[p_down->first], best_cost);
                if (cost < best_cost)
                {
                    best_cost = cost;
                    best_p = p_down;

                    if (best_cost == 0)
                        break;
                }

                --p_down;
                dif_down = abs(p->second - p_down->second);
            }
        }

        perm.push_back(best_p->first);
        density_list.erase(p);
        p = best_p;

        cost_perm += best_cost;
        --n_cur_samples;
    }	
	
    uint64_t cost_tsp2 = 0;
    for (int i = 1; i < n_h_samples; ++i)
        cost_tsp2 += bit_cost(mc_vectors[perm[i - 1]], mc_vectors[perm[i]]);
    
    // Perform permutation
    int64_t no1_orig = 0;
    int64_t no1_perm = 0;
    for (int64_t i = id_start*s->vec_len; i < (int64) (id_stop * s->vec_len); ++i)
        for (int x = data[i]; x; x &= x - 1)
            ++no1_orig;
    
    // Prepare permutation vector
    v_perm.clear();
    v_perm.resize(n_h_samples, 0);
    for (int i = 0; i < n_h_samples; ++i)
        v_perm[perm[i]] = i;
    
    uchar_t *old_vec = new uchar_t[s->vec_len];
    
    for (uint64_t i = id_start; i < id_stop; ++i)
    {
        auto new_vec = data + i*s->vec_len;
        memcpy(old_vec, new_vec, s->vec_len);
        
        fill_n(new_vec, s->vec_len, 0);
        
        uint32 x;
        for (x = 0; x + 8 < (uint32) n_h_samples;)
        {
            int x8 = x / 8;
            int d_x8 = old_vec[x8];
            if(!d_x8)
            {
                x += 8;
                continue;
            }
		
            int x_p = 1 << 7;
            for(int i = 0; i < 8; ++i, ++x, x_p >>= 1)
            {
                auto j = v_perm[x];
                uchar_t j_p = perm_lut8[j % 8];
                int j8 = j / 8;
                
                if (d_x8 & x_p)
                    new_vec[j8] += j_p;
            }
        }
		
        int x8 = x / 8;
        int d_x8 = old_vec[x8];

        for (; x < (uint32) n_h_samples; ++x)
        {
            auto j = v_perm[x];
            uchar_t x_p = perm_lut8[x % 8];
            uchar_t j_p = perm_lut8[j % 8];
            int j8 = j / 8;

            if (d_x8 & x_p)
                new_vec[j8] += j_p;
        }
    }
    
    delete[] old_vec;
    
    // Calculate no. of 1s to verify the correctness of permutation - to be removed in the final version !!!
    for (int64_t i = id_start*(s->vec_len); i < (int64) (id_stop * s->vec_len); ++i)
        for (int x = data[i]; x; x &= x - 1)
            ++no1_perm;
    
    if (no1_orig != no1_perm)
        cout << "\nWrong permutation!  " << no1_orig << " " << no1_perm << endl;
}
