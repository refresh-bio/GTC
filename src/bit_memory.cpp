/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include "defs.h"
#include "bit_memory.h"
#include <algorithm>

using namespace std;

// ********************************************************************************************
CBitMemory::CBitMemory()
{
    mem_buffer	     = NULL;
    mem_buffer_size  = 0;
    mem_buffer_pos   = 0;
    word_buffer_pos  = 0;
    word_buffer		 = 0;
    word_buffer_size = 32;
    mode			 = mode_none;
    mem_buffer_ownership = false;
    
    for(int32 i = 0; i < 32; ++i)
        n_bit_mask[i] = (1u << i) - 1;
}

// ********************************************************************************************
CBitMemory::CBitMemory(const CBitMemory &y)
{
    CBitMemory &x = const_cast<CBitMemory &>(y);
    
    mem_buffer = x.mem_buffer;
    x.mem_buffer = NULL;
    mem_buffer_ownership = x.mem_buffer_ownership;
    x.mem_buffer_ownership = false;
    
    mem_buffer_pos = x.mem_buffer_pos;
    x.mem_buffer_pos = 0;
    
    mem_buffer_size = x.mem_buffer_size;
    x.mem_buffer_size = 0;
    
    word_buffer = x.word_buffer;
    word_buffer = 0;
    
    word_buffer_pos = x.word_buffer_pos;
    x.word_buffer_pos = 0;
    
    word_buffer_size = x.word_buffer_size;
    
    mode = x.mode;
    mode = mode_none;
    
    copy_n(x.n_bit_mask, 32, n_bit_mask);
}

// ********************************************************************************************
CBitMemory::~CBitMemory()
{
    if(mem_buffer && mem_buffer_ownership)
        delete[] mem_buffer;
}

// ********************************************************************************************
bool CBitMemory::TakeOwnership()
{
    if(mem_buffer && mem_buffer_ownership)
    {
        mem_buffer_ownership = false;
        return true;
    }
    
    return false;
}

// ********************************************************************************************
bool CBitMemory::Open(uchar *p, int64 size, bool force_open)
{
    if(!force_open)
    {
        if(mode != mode_none)
            return false;
    }
    
    if(mem_buffer && mem_buffer_ownership)
        delete[] mem_buffer;
    
    mem_buffer_size = size;
    /*	if(!mem_buffer_size)
     mem_buffer_size = 1;
     mem_buffer = new uchar[mem_buffer_size];
     copy_n(p, size, mem_buffer);*/
    mem_buffer = p;
    mem_buffer_ownership = false;

    if(!mem_buffer_size)
    {
        mem_buffer_size = 1;
        mem_buffer = new uchar[mem_buffer_size];
        mem_buffer_ownership = true;
        copy_n(p, size, mem_buffer);
    }
    
    mode = mode_mem_read;
    
    word_buffer_size = 8;
    mem_buffer_pos   = 0;
    
    return mode == mode_mem_read;
}

// ********************************************************************************************
bool CBitMemory::Restart()
{
    mode = mode_mem_read;
    
    word_buffer_size = 32;
    mem_buffer_pos   = 0;
    word_buffer_pos  = 0;
    word_buffer		 = 0;
    
    return true;
}

// ********************************************************************************************



// ********************************************************************************************
bool CBitMemory::Create(int64 size)
{
    if(mode != mode_none)
        return false;
    
    if(mem_buffer && mem_buffer_ownership)
        delete[] mem_buffer;
    
    if(!size)
        size = 1;
    mem_buffer_size = size;
    mem_buffer = new uchar[mem_buffer_size];
    mem_buffer_ownership = true;
    
    mode = mode_mem_write;
    mem_buffer_pos = 0;
    
    word_buffer_size = 32;
    
    return mode == mode_mem_write;
}

// ********************************************************************************************
bool CBitMemory::Close()
{
    if(mode != mode_mem_write && mode != mode_mem_read)
        return false;
    
    if(mode == mode_mem_write)
    {
        if(word_buffer_pos)
            FlushPartialWordBuffer();
    }
    
    if(mem_buffer && mem_buffer_ownership)
        delete[] mem_buffer;
    
    mem_buffer = NULL;
    mem_buffer_ownership = false;
    mode = mode_none;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::Complete()
{
    if(mode != mode_mem_write && mode != mode_mem_read)
        return false;
    
    if(mode == mode_mem_write)
    {
        if(word_buffer_pos)
            FlushPartialWordBuffer();
    }
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::SetPos(int64 pos)
{
    if(mode != mode_file_read && mode != mode_mem_read)
        return false;
    
    if((int64) pos > mem_buffer_size)
        return false;
    mem_buffer_pos = pos;
    
    word_buffer_pos = 0;
    word_buffer     = 0;
    
    return true;
}
