/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include <algorithm>
#include "huffman.h"
#include "defs.h"

using namespace std;

// ********************************************************************************************
uint32 int_log(uint32 x, uint32 base)
{
    uint32 r = 0;
    
    if(base == 0)
        return 1;
    if(base == 1)
        base++;
    
    for(uint32 tmp = base; tmp <= x; tmp *= base)
        ++r;
    
    return r;
}

// ********************************************************************************************
CHuffman::CHuffman(uint32 _size)
{
    size = _size;
    if(size)
    {
        tree  = new t_node[2*size-1];
        codes = new t_code[2*size-1];
        heap  = new t_frequency[size];
    }
    else
    {
        tree  = NULL;
        codes = NULL;
        heap  = NULL;
    }
    n_symbols = 0;
    min_len   = 1;
    
    speedup_tree = NULL;
    bit_memory = NULL;
}

// ********************************************************************************************
CHuffman::~CHuffman()
{
    if(tree)
        delete[] tree;
    if(codes)
        delete[] codes;
    if(heap)
        delete[] heap;
    
    if(speedup_tree)
        delete[] speedup_tree;
    
    if(bit_memory)
        delete bit_memory;
    
    if(speedup_lut)
        delete[] speedup_lut;
}

// ********************************************************************************************
bool CHuffman::Restart(uint32 _size)
{
    if(tree)
        delete[] tree;
    if(codes)
        delete[] codes;
    if(heap)
        delete[] heap;
    if(speedup_tree)
        delete[] speedup_tree;
    if(speedup_lut)
        delete[] speedup_lut;
    
    size = _size;
    if(size)
    {
        tree  = new t_node[2*size-1];
        codes = new t_code[2*size-1];
        heap  = new t_frequency[size];
    }
    else
    {
        tree  = NULL;
        codes = NULL;
        heap  = NULL;
    }
    n_symbols = 0;
    
    speedup_tree = NULL;
    
    return true;
}

// ********************************************************************************************
bool CHuffman::RestartDecompress(uint32 _size, uint32 _root_id)
{
    if(tree)
        delete[] tree;
    if(codes)
        delete[] codes;
    if(heap)
        delete[] heap;
    if(speedup_tree)
        delete[] speedup_tree;
    
    size = _size;
    if(size)
    {
        tree  = new t_node[_root_id - _size + 2];
        codes = NULL;
        heap = NULL;
    }
    else
    {
        tree  = NULL;
        codes = NULL;
        heap  = NULL;
    }
    n_symbols = size;
    
    speedup_tree = NULL;
    
    return true;
}

// ********************************************************************************************
CHuffman::t_code* CHuffman::Complete(bool compact)
{
    int32 i;
    
    if(!n_symbols)
        return NULL;
    
    // Make heap of symbols
    make_heap(heap, heap+n_symbols);
    
    // Prepare leaves of the tree
    for(i = 0; i < n_symbols; ++i)
    {
        codes[i].code = 0;
        codes[i].len  = 0;
        tree[i].left_child = -1;
        tree[i].right_child = -1;
    }
    for(i = n_symbols; i < 2*n_symbols-1; ++i)
    {
        codes[i].code = 0;
        codes[i].len  = 0;
    }
    
    // Build tree
    int32 heap_size = n_symbols;
    // Remove symbols with 0 frequency
    if(compact)
        while(heap_size > 2 && heap[0].frequency == 0)
            pop_heap(heap, heap+heap_size--);
    
    int32 present_symbols = heap_size;
    
    if(!present_symbols)
        return codes;
    
    for(i = 0; i < present_symbols-1; ++i)
    {
        t_frequency left = heap[0];
        pop_heap(heap, heap+heap_size--);
        t_frequency right = heap[0];
        pop_heap(heap, heap+heap_size--);
        
        heap[heap_size].symbol = n_symbols+i;
        heap[heap_size].frequency = left.frequency + right.frequency;
        push_heap(heap, heap+ ++heap_size);
        
        tree[n_symbols+i].left_child  = left.symbol;
        tree[n_symbols+i].right_child = right.symbol;
    }
    
    // Compute codes
    for(i = n_symbols+present_symbols-2; i >= n_symbols; --i)
    {
        codes[tree[i].left_child].len   = codes[i].len+1;
        codes[tree[i].left_child].code  = (codes[i].code << 1);
        codes[tree[i].right_child].len  = codes[i].len+1;
        codes[tree[i].right_child].code = (codes[i].code << 1) | 1;
    }
    
    root_id = n_symbols + present_symbols - 2;
    cur_id = root_id;
    
    return codes;
}

// ********************************************************************************************
bool CHuffman::StoreTree(uchar *&mem, uint32 &len)
{
    if(bit_memory)
        delete bit_memory;
    
    bit_memory = new CBitMemory();
    
    bits_per_id = int_log(size, 2);
    if(size & (size-1))			// size is not power of 2
        bits_per_id++;
    
    min_len = 32;
    for(int i = 0; i < n_symbols; ++i)
        if(codes[i].len < min_len && codes[i].len > 0)
            min_len = codes[i].len;
    
    max_len = 0;
    for(int i = 0; i < n_symbols; ++i)
        if(codes[i].len > max_len)
            max_len = codes[i].len;
    
    if(n_symbols == 1 || min_len > max_len)
        min_len = 0;
    
    int32 node_id = root_id;
    bit_memory->PutWord(root_id);
    bit_memory->PutWord(n_symbols);
    bit_memory->PutByte((uchar) min_len);
    bit_memory->PutByte((uchar) max_len);
    EncodeProcess(node_id);
    bit_memory->FlushPartialWordBuffer();
    
    mem = bit_memory->mem_buffer;
    len = (int32) bit_memory->mem_buffer_pos;
    
    return true;
}

// ********************************************************************************************
void CHuffman::EncodeProcess(int32 node_id)
{
    if(tree[node_id].left_child == -1)		// leaf
    {
        bit_memory->PutBit(1);
        bit_memory->PutBits(node_id, bits_per_id);
    }
    else
    {
        bit_memory->PutBit(0);
        EncodeProcess(tree[node_id].left_child);
        EncodeProcess(tree[node_id].right_child);
    }
}

// ********************************************************************************************
bool CHuffman::LoadTree(uchar *mem, uint32 len)
{
    if(bit_memory)
        delete bit_memory;
    
    bit_memory = new CBitMemory();
    bit_memory->Open(mem, len);
    
    uint32 tmp;
    
    bit_memory->GetWord(tmp);
    root_id = tmp;
    bit_memory->GetWord(tmp);
    n_symbols = tmp;
    tmp_id = root_id;
    cur_id = root_id;
    
    bit_memory->GetByte(min_len);
    bit_memory->GetByte(max_len);
    
    RestartDecompress(n_symbols, root_id);
    bits_per_id = int_log(size, 2);
    if(size & (size-1))			// size is not power of 2
        bits_per_id++;
    
    
    int32 node_id = root_id - n_symbols + 1;
    root_id = node_id;
    tmp_id = root_id;
    cur_id = root_id;
    DecodeProcess(node_id);
    
    if(!min_len)
        min_len = 1;
    
    ComputeSpeedupTree();
    ComputeSpeedupTreeLongest();
    
    if(min_len > max_len)
        min_len = 0;
    
    delete bit_memory;
    bit_memory = NULL;
    
    return true;
}

// ********************************************************************************************
int32 CHuffman::DecodeProcess(int32 node_id)
{
    uint32 flag = 0;
    uint32 tmp;
    
    bit_memory->GetBit(flag);
    
    if(!flag)
    {
        --tmp_id;
        tree[node_id].left_child  = DecodeProcess(tmp_id);
        tree[node_id].right_child = DecodeProcess(tmp_id);
        return node_id;
    }
    else
    {
        bit_memory->GetBits(tmp, bits_per_id);
        return -tmp;
    }
    
}

// ********************************************************************************************
void CHuffman::ComputeSpeedupTree()
{
    int i, j;
    
    if(!min_len)
        return;
    
    if(speedup_tree)
        delete[] speedup_tree;
    speedup_tree = new int32[(uint32) (1 << min_len)];
    
    for(i = 0; i < (1 << min_len); ++i)
    {
        cur_id = root_id;
        for(j = min_len-1; j >= 0; --j)
            Decode(i & (1 << j));
        
        speedup_tree[i] = cur_id;
    }
    
    cur_id = root_id;
    tmp_id = root_id;
}

// ********************************************************************************************
void CHuffman::ComputeSpeedupTreeLongest()
{
    int i, j;
    
    if(!max_len)
        return;
    
    if(speedup_lut)
        delete[] speedup_lut;
    
    uint32_t size;
    max_lut_len = MAX_HUF_LUT_LEN;
    
    if(max_len <= max_lut_len)
    {   size = (uint32_t) (1 << max_len);
        max_lut_len = max_len;
    }
    else
        size = (uint32_t) (1 << max_lut_len);
    
    speedup_lut = new uint16_t[size];
    
    for(i = 0; i < (int) size; ++i)
    {
        cur_id = root_id;
        
        for(j = max_lut_len-1; j >= ((int)max_lut_len-(int)min_len); --j)
            Decode(i & (1 << j));
        
        
        while (cur_id > 0 && j >=0)
        {
            Decode(i & (1 << j));
            --j;
        }
        char code_len = max_lut_len-j-1;
        
        uint32_t end_pos = i;
        for(int idx = j; idx >= 0; idx--)
            end_pos |= (1<<idx);
        
        uint16_t tuple;
        // [ 10 bits      | 1 bit                                  | 5bits          ]
        // [ abs(cur_id) | flag=1, if positive curr_id (not found) | number_of_bits ]
        if(cur_id<=0)
        {
            tuple = ((-cur_id << 6)|(uint8_t)code_len);
        }
        else
        {
            tuple = ((cur_id << 6)|(uint8_t)code_len);
            tuple = tuple | (1 << 5);
        }
        
        for(int idx = i; idx <= (int) end_pos; idx++)
        {
            speedup_lut[idx] = tuple;
        }
        i = end_pos;
    }
    
    cur_id = root_id;
    tmp_id = root_id;
}
