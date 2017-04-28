/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#ifndef _HUFFMAN_H
#define _HUFFMAN_H

#include "defs.h"
#include "bit_memory.h"

class CHuffman {
    int32 size;
    
    typedef struct tag_frequency {
        uint32 symbol;
        uint32 frequency;
        
        bool operator<(const struct tag_frequency &x) const {
            return frequency > x.frequency || (frequency == x.frequency && symbol > x.symbol);
        }
    } t_frequency;
    
public:
    typedef struct {
        uint32 code;
        uint32 len;
    } t_code;
    
private:
public:
    typedef struct {
        int32 left_child = 0;
        int32 right_child = 0;
    } t_node;
    
    t_frequency *heap  = NULL;
    t_node *tree  = NULL;
    int32 root_id;
    int32 n_symbols;
    int32 cur_id;
    int32 tmp_id;
    uint32 min_len, max_len, max_lut_len = 0;
    
    CBitMemory *bit_memory  = NULL;
    uint32 bits_per_id;
    
    int32 *speedup_tree = NULL;
    uint16_t *speedup_lut = NULL;
    
    void EncodeProcess(int32 node_id);
    int32 DecodeProcess(int32 node_id);
    void ComputeSpeedupTree();
    void ComputeSpeedupTreeLongest();
    
public:
    t_code *codes;
    
    CHuffman(uint32 _size = 0);
    ~CHuffman();
    
    bool Restart(uint32 _size = 0);
    bool RestartDecompress(uint32 _size = 0, uint32 _root_id = 0);
    inline bool Insert(const uint32 frequency);
    t_code* Complete(bool compact = true);
    inline int32 Decode(const uint32 bit);
    inline int32 DecodeFast(const uint32 bits);
    inline int32 DecodeFastLut(const uint32 bits, uint8_t & no_bits);
    
    bool StoreTree(uchar *&mem, uint32 &len);
    bool LoadTree(uchar *mem, uint32 len);
};

// ********************************************************************************************
bool CHuffman::Insert(const uint32 frequency)
{
    if(n_symbols == size)
        return false;
    
    heap[n_symbols].symbol    = n_symbols;
    heap[n_symbols].frequency = frequency;
    n_symbols++;
    
    return true;
}

// ********************************************************************************************
int32 CHuffman::Decode(const uint32 bit)
{
    if(cur_id <= 0)
        cur_id = root_id;
    if(bit)
        cur_id = tree[cur_id].right_child;
    else
        cur_id = tree[cur_id].left_child;

    if(cur_id <= 0)
        return -cur_id;				// Symbol found
    else
        return -1;					// Not found yet
}

// ********************************************************************************************
inline int32 CHuffman::DecodeFast(const uint32 bits)
{
    cur_id = speedup_tree[bits];

    if(cur_id <= 0)
        return -cur_id;				// Symbol found
    else
        return -1;					// Not found yet
}

// ********************************************************************************************
inline int32 CHuffman::DecodeFastLut(const uint32 bits, uint8_t & no_bits)
{
    uint16_t tuple = speedup_lut[bits];
    // [ 10 bits      | 1 bit                             | 5bits          ]
    // [ abs(cur_id) | flag = 1, if positive (not found) | number_of_bits ]
    
    no_bits = tuple & 31;
    
    cur_id = (tuple >> 6);
    if((tuple & 32) > 0) // positive, not found
    {
        return -1;
    }
    else
    {
        return cur_id;
    }
}

#endif
