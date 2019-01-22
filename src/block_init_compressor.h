/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */
#ifndef block_init_compressor_h
#define block_init_compressor_h

#include <iostream>
#include "params.h"
#include "defs.h"
#include "bit_memory.h"
#include "compression_settings.h"
#include <array>
#include "nmmintrin.h"
#include <random>
#include <cstring>

typedef array<uint64_t, MC_ARRAY_SIZE> mc_vec_t;

// Count number of bits between two arrays that are different (same positions)
// Operation POPCNT is used (available for modern CPUs)
template <typename T>
inline uint64_t bit_cost(const T &x, const T&y)
{
    uint64_t r = 0;

    int max_i = x.size();

    switch (max_i % 4)
    {
    case 3:
        --max_i;
        r += _mm_popcnt_u64(x[max_i] ^ y[max_i]);
    case 2:
        --max_i;
        r += _mm_popcnt_u64(x[max_i] ^ y[max_i]);
    case 1:
        --max_i;
        r += _mm_popcnt_u64(x[max_i] ^ y[max_i]);
    }

    for (int i = max_i - 1; i >= 0;)
    {
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
    }

    return r;
}

template <typename T>
inline uint64_t bit_cost(const T &x, const T&y, uint64_t best_cost)
{
    uint64_t r = 0;
    int max_i = x.size();

    switch (max_i % 4)
    {
    case 3:
        --max_i;
        r += _mm_popcnt_u64(x[max_i] ^ y[max_i]);
    case 2:
        --max_i;
        r += _mm_popcnt_u64(x[max_i] ^ y[max_i]);
    case 1:
        --max_i;
        r += _mm_popcnt_u64(x[max_i] ^ y[max_i]);
    }

    for (int i = max_i - 1; i >= 0 && r < best_cost;)
    {
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
    }

    return r;
}

class BlockInitCompressor{

    CBitMemory bm;
    
    CompSettings * s = nullptr;

    uint64_t cur_no_vec;
    
    uint64_t comp_no_matches, comp_no_literals, comp_zero_run, comp_ones_run, same_vec_match, n_copies;
     uint64_t cur_vec_id;
    uint32_t uniq_vec_counter;
    
    uint32_t no_non_copy; //no_vec - rank_zeros_only_vector(no_vec) - rank_copy_vector(no_vec);
    uint32_t no_copy; //rank_copy_vector(no_vec);
    uint32_t no_zeros_only; //rank_zeros_only_vector(no_vec);
   
    uchar_t * data = nullptr;
    uchar_t * counters = nullptr;
    uchar_t * temp_vec = nullptr;
    
    uint32_t * comp_pos_non_copy = nullptr;
    uint32_t * comp_pos_copy = nullptr;
    
    uint32_t perm_lut8[8];
    uint64_t perm_lut64[64];
    
    //hashers
    uint64_t  ht_checks;
    uint32_t  **ht = nullptr;
    uint32_t  *ht_size = nullptr;
    uint32_t  *ht_fill = nullptr;
    uint32_t  *ht_entries = nullptr;
    uint32_t *vec_hash = nullptr;
    uint64_t vec_hash_size;
    
    vector<bool> * zeros_only = nullptr;
    vector<bool> * copy = nullptr;
    
    // Hashing functions; it calculates hash for position pointed by in_vec_pos
    uint64_t hash_fun(uint64_t vec_id, uint64_t in_vec_pos, uint64_t len, uint64_t size);
    
    // Insert into HT
    bool ht_insert(uint64_t vec_id, uint64_t in_vec_pos);
    
    // Remove from HT
    bool ht_remove(uint64_t vec_id, uint64_t in_vec_pos);
    
    // Resize HT
    void restruct(uint64_t in_vec_pos);
    
    // Inserting to hash table for whole vectors
    // Resizing is not needed  - initial size is more than enough
    bool ht_insert_vec(uint64_t vec_id, uint32_t key);

    uint32_t idx_oldest_vec_in_ht_parts = 0;
    uint32_t idx_oldest_vec_in_ht_vecs  = 0;
    
    uchar_t lookup_t_ones[256] = {0};
    
    void Allocate();
    void Deallocate();
    bool allocated;
    void Clear();
    void comp_vec(uint64_t vec_id);
    void encode_literal_run(int32_t len, const uchar * lit_run);
    char bits_used(unsigned int n);
    uchar_t get_ones_group(uint64_t vec_id);
    void permute_range_vec(uint64_t id_start, uint64_t id_stop, vector<int> &v_perm);
    
    double aux_dividor; //help for get_ones_group 
    
    uint32_t n_vec_in_ht_parts;
    uint32_t n_vec_in_ht_vecs;
    
public:
    
    BlockInitCompressor(CompSettings * _settings)
    {
        s = _settings;
        allocated = false;
        aux_dividor = (double)s->vec_len*BITS_IN_BYTE/s->ones_ranges;

        no_non_copy = 0;
        no_copy = 0;
        no_zeros_only = 0;
        uniq_vec_counter = 0;
        n_vec_in_ht_parts = 0;
        n_vec_in_ht_vecs  = 0;
        
        // Permutation LUT
        for(int i = 0; i < 8; ++i)
            perm_lut8[i] = 1 << (7 - i);	

        for(int i = 0; i < 64; ++i)
            perm_lut64[i] = 1ull << i;	
    }
    ~BlockInitCompressor()
    {
        Deallocate();
    }
    
    void SetBlock(uint64_t _cur_no_vec, uchar_t * _data);
    void PermuteBlock(vector<int> & perm, bool permute = true);
    bool Compress(vector<bool> &zeros, vector<bool> &copies, uchar * &compressedBlock, size_t & compressed_size, uint32_t *& origin_of_copy);
};

#endif /* block_init_compressor_h */


