/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */
#ifndef _BIT_MEMORY_H
#define _BIT_MEMORY_H

#include "defs.h"
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

typedef enum {mode_none, mode_file_read, mode_file_write, mode_file_read_ra, mode_mem_read, mode_mem_write, mode_mem_read_ra} t_mode;

// ********************************************************************************************
class CBitMemory {
public:
    uchar *mem_buffer;
    bool mem_buffer_ownership;
    int64 mem_buffer_pos;
    int32 word_buffer_pos;
private:
    int64 mem_buffer_size;
    uint32 word_buffer;
    int32 word_buffer_size;
    t_mode mode;
    
public:
    uint32 n_bit_mask[32];
    
    inline bool FlushFullWordBuffer();
    inline bool FlushPartialWordBuffer();
    inline bool FlushInputWordBuffer();
    
    CBitMemory();
    CBitMemory(const CBitMemory &y);
    ~CBitMemory();
    
    bool Open(uchar *p, int64 size, bool force_open = false);
    bool Create(int64 size = 1);
    bool Complete();
    bool Close();
    bool Restart();
    bool RestartWrite();
    bool TakeOwnership();
    
    inline uint64 GetPos(void);
    bool SetPos(int64 pos);
    
    inline int32 GetWordPos(void);
    
    inline bool PutBit(const uint32 word);
    inline bool Put2Bits(const uint32 word);
    inline bool PutBits(uint32 word, int32 n_bits);
    inline bool PutBytes(const unsigned char *data, int64 n_bytes);
    inline bool PutByte(const unsigned char byte);
    inline bool PutNBytes(const unsigned char byte, int64 n_bytes);
    inline bool PutWord(const uint32 data);
    inline bool PutDWord(const uint64 data);
    inline bool Put2Bytes(const uint32 data);
    
    inline bool GetBit(uint32 &word);
    inline bool Get2Bits(uint32 &word);
    inline bool GetBits(uint32 &word, uint32 n_bits);
    inline bool GetBytes(unsigned char *data, uint64 n_bytes);
    inline bool GetBool(bool &byte);
    inline bool GetByte(uint32 &byte);
    inline bool GetWord(uint32 &data);
    inline bool GetWord(int32 &data);
    inline bool GetDWord(uint64 &data);
    inline bool Get2Bytes(uint32 &data);
    
    inline bool SkipBytes(uint32 n_bytes);
    
    inline uint32 BitLength(const uint64 x);
    
    uint64 GetFilePos() {return mem_buffer_pos;};
    
    inline bool UndoGetByte(uint32 &byte);//AD //while getting
    inline bool UndoGetBits(uint32 n_bits); //AD //while getting
    inline bool discardBits(uint32 word);  //AD //while putting
    inline bool GetBitsAndDiscard(uint32 n_bits); //while getting

};


// ********************************************************************************************
bool CBitMemory::GetBitsAndDiscard(uint32 n_bits)
{
    uint32 count;
   
    while(n_bits)
    {
        if(word_buffer_pos == 0)
        {
            count = n_bits/8;
            if(count > 1)
            {
                if(mem_buffer_pos + count - 1 >= mem_buffer_size)
                    return false;
                //  byte = mem_buffer[mem_buffer_pos++];
                mem_buffer_pos += count - 1;
                n_bits -= (count - 1) << 3;  //*8
            }
            if(!GetByte(word_buffer))
                return false;
            word_buffer_pos = 8;
        }
        
        if((int32) n_bits > word_buffer_pos)
        {
          //  word <<= word_buffer_pos;
          //  word += word_buffer & n_bit_mask[word_buffer_pos];
            n_bits -= word_buffer_pos;
            word_buffer_pos = 0;
        }
        else
        {
          //  word <<= n_bits;
            word_buffer_pos -= n_bits;
          //  word += (word_buffer >> word_buffer_pos) & n_bit_mask[n_bits];
            return true;
        }
    }
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::discardBits(uint32 n_bits)  //while putting
{
    uint32_t word;
    while(n_bits)
    {
        if((int32) n_bits <= word_buffer_pos)
        {
            word_buffer_pos -= n_bits;
            word_buffer = word_buffer >> n_bits ;
            if(!word_buffer_pos)
                word_buffer = 0;
            return true;
        }
        else //(n_bits > word_buf_pos)
        {
            n_bits -= word_buffer_pos;
            word_buffer_pos = 0;
            word_buffer = 0;
        }
        
        mem_buffer_pos = mem_buffer_pos - 4;
        //memcpy(&word_buffer, &mem_buffer[mem_buffer_pos], 4);
        GetWord(word);
        mem_buffer_pos = mem_buffer_pos - 4;
        word_buffer = word;
        word_buffer_pos = 32;

    }
    return false;
}


// ********************************************************************************************
uint64 CBitMemory::GetPos(void)
{
    return mem_buffer_pos;
}

// ********************************************************************************************
int32 CBitMemory::GetWordPos(void)
{
    return word_buffer_pos;
}


// ********************************************************************************************
bool CBitMemory::PutBit(const uint32 word)
{
    if(word_buffer_pos < word_buffer_size)
    {
        word_buffer <<= 1;
        word_buffer += word;
        ++word_buffer_pos;
    }
    else
    {
        PutWord(word_buffer);
        word_buffer_pos = 1;
        word_buffer = word;
    }
    
    return true;
};

// ********************************************************************************************
bool CBitMemory::Put2Bits(const uint32 word)
{
    if(word_buffer_pos + 2 <= word_buffer_size)
    {
        word_buffer <<= 2;
        word_buffer += word;
        word_buffer_pos += 2;
    }
    else if(word_buffer_pos == word_buffer_size)
    {
        PutWord(word_buffer);
        word_buffer_pos = 2;
        word_buffer = word;
    }
    else
    {
        word_buffer <<= 1;
        word_buffer += word >> 1;
        PutWord(word_buffer);
        word_buffer = word & 1;
        word_buffer_pos = 1;
    }
    
    return true;
};

// ********************************************************************************************
bool CBitMemory::PutBits(uint32 word, int32 n_bits)
{
    int32 rest_bits = word_buffer_size - word_buffer_pos;
    if(n_bits >= rest_bits)
    {
        n_bits -= rest_bits;
        word_buffer <<= rest_bits;
        word_buffer += word >> n_bits;
        word &= n_bit_mask[n_bits];
        word_buffer_pos = 0;
        PutWord(word_buffer);
        word_buffer = 0;
    }
    
    word_buffer     <<= n_bits;
    word_buffer     += word;
    word_buffer_pos += n_bits;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::FlushFullWordBuffer()
{
    PutWord(word_buffer);
    
    word_buffer     = 0;
    word_buffer_pos = 0;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::FlushPartialWordBuffer()
{
    word_buffer <<= (32 - word_buffer_pos) & 7;
    
    if(word_buffer_pos > 24)
        PutByte(word_buffer >> 24);
    if(word_buffer_pos > 16)
        PutByte((word_buffer >> 16) & 0xFF);
    if(word_buffer_pos > 8)
        PutByte((word_buffer >> 8) & 0xFF);
    if(word_buffer_pos > 0)
        PutByte(word_buffer & 0xFF);
    
    word_buffer     = 0;
    word_buffer_pos = 0;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::FlushInputWordBuffer()
{
    word_buffer_pos = 0;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::PutByte(const unsigned char byte)
{
    if(mem_buffer_pos + 1 > mem_buffer_size)
    {
        mem_buffer_size = (uint64) ((mem_buffer_pos + 1) * 1.5);
        uchar *new_mem_buffer = new uchar[mem_buffer_size];
        if(mem_buffer)
        {
            copy_n(mem_buffer, mem_buffer_pos, new_mem_buffer);
            delete[] mem_buffer;
        }
        mem_buffer = new_mem_buffer;
    }
    
    mem_buffer[mem_buffer_pos++] = byte;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::Put2Bytes(const uint32 data)
{
    PutByte((unsigned char) (data >> 8));
    PutByte((unsigned char) (data & 0xFF));
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::PutBytes(const unsigned char *data, int64 n_bytes)
{
   
    if(mem_buffer_pos + n_bytes > mem_buffer_size)
    {
        mem_buffer_size = (uint64) ((mem_buffer_pos + n_bytes) * 1.5);
        uchar *new_mem_buffer = new uchar[mem_buffer_size];
        if(mem_buffer)
        {
            copy_n(mem_buffer, mem_buffer_pos, new_mem_buffer);
            delete[] mem_buffer;
        }
        mem_buffer = new_mem_buffer;
    }
    copy_n(data, n_bytes, mem_buffer+mem_buffer_pos);
    mem_buffer_pos += n_bytes;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::PutNBytes(const unsigned char byte, int64 n_bytes)
{
    if(mem_buffer_pos + n_bytes > mem_buffer_size)
    {
        mem_buffer_size = (uint64) ((mem_buffer_pos + n_bytes) * 1.5);
        uchar *new_mem_buffer = new uchar[mem_buffer_size];
        if(mem_buffer)
        {
            copy_n(mem_buffer, mem_buffer_pos, new_mem_buffer);
            delete[] mem_buffer;
        }
        mem_buffer = new_mem_buffer;
    }
    fill_n(mem_buffer+mem_buffer_pos, n_bytes, byte);
    mem_buffer_pos += n_bytes;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::PutDWord(const uint64 data)
{
    PutByte(data >> 56);
    PutByte((data >> 48) & 0xFF);
    PutByte((data >> 40) & 0xFF);
    PutByte((data >> 32) & 0xFF);
    PutByte((data >> 24) & 0xFF);
    PutByte((data >> 16) & 0xFF);
    PutByte((data >> 8) & 0xFF);
    PutByte(data & 0xFF);
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::PutWord(const uint32 data)
{
    PutByte(data >> 24);
    PutByte((data >> 16) & 0xFF);
    PutByte((data >> 8) & 0xFF);
    PutByte(data & 0xFF);
    
    return true;
}
// ********************************************************************************************
uint32 CBitMemory::BitLength(const uint64 x)
{
    for(uint32 i = 0; i < 32; ++i)
        if(x < (1ull << i))
            return i;
    
    return 64;
}

// ********************************************************************************************
bool CBitMemory::GetBit(uint32 &word)
{
    if(word_buffer_pos == 0)
    {
        if(!GetByte(word_buffer))
            return false;
        word_buffer_pos = 7;
        word = word_buffer >> 7;
    }
    else
        word = (word_buffer >> (--word_buffer_pos)) & 1;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::Get2Bits(uint32 &word)
{
    if(word_buffer_pos >= 2)
    {
        word = (word_buffer >> (word_buffer_pos-2)) & 3;
        word_buffer_pos -= 2;
    }
    else if(word_buffer_pos == 0)
    {
        if(!GetByte(word_buffer))
            return false;
        word = word_buffer >> 6;
        word_buffer_pos = 6;
    }
    else
    {
        word = (word_buffer & 1) << 1;
        if(!GetByte(word_buffer))
            return false;
        word += word_buffer >> 7;
        word_buffer_pos = 7;
    }
    
    return true;
};

// ********************************************************************************************
bool CBitMemory::GetBits(uint32 &word, uint32 n_bits)
{
    word = 0;
    while(n_bits)
    {
        if(word_buffer_pos == 0)
        {
            if(!GetByte(word_buffer))
                return false;
            word_buffer_pos = 8;
        }
        
        if((int32) n_bits > word_buffer_pos)
        {
            word <<= word_buffer_pos;
            word += word_buffer & n_bit_mask[word_buffer_pos];
            n_bits -= word_buffer_pos;
            word_buffer_pos = 0;
        }
        else
        {
            word <<= n_bits;
            word_buffer_pos -= n_bits;
            word += (word_buffer >> word_buffer_pos) & n_bit_mask[n_bits];
            return true;
        }
    }
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::UndoGetBits(uint32 n_bits)
{
    
    while(n_bits)
    {
        if(word_buffer_pos == 8)
        {
            if(!UndoGetByte(word_buffer))
                return false;
            word_buffer_pos = 0;
        }
        
        
        if((int32) n_bits <= 8 - word_buffer_pos)
        {
            word_buffer_pos += n_bits;
            return true;
        }
        else
        {
            n_bits -= (8 - word_buffer_pos);
            word_buffer_pos = 8;
        }
        
        
    }
    
    return true;
}



// ********************************************************************************************
bool CBitMemory::GetByte(uint32 &byte)
{
    if(mem_buffer_pos >= mem_buffer_size)
        return false;
    byte = mem_buffer[mem_buffer_pos++];
    return true;
};



// ********************************************************************************************
bool CBitMemory::UndoGetByte(uint32 &byte)
{
    if(mem_buffer_pos <= 1)
        return false;
    --mem_buffer_pos;
    --mem_buffer_pos;
    byte = mem_buffer[mem_buffer_pos++];
    
    return true;
};

// ********************************************************************************************
bool CBitMemory::GetBool(bool &byte)
{
    if(mem_buffer_pos >= mem_buffer_size)
        return false;
    byte = (mem_buffer[mem_buffer_pos++] == 1);
    
    return true;
};

// ********************************************************************************************
bool CBitMemory::GetBytes(unsigned char *data, uint64 n_bytes)
{
    if(mem_buffer_pos + (int64) n_bytes >= mem_buffer_size)
        return false;
    
    copy_n(mem_buffer+mem_buffer_pos, n_bytes, data);
    mem_buffer_pos += n_bytes;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::GetWord(uint32 &data)
{
    uint32 c = 0;
    bool r;
    
    r = GetByte(c);
    data = c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    
    return r;
}

// ********************************************************************************************
bool CBitMemory::GetWord(int32 &data)
{
    uint32 c;
    bool r;
    
    r = GetByte(c);
    data = c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    
    return r;
}

// ********************************************************************************************
bool CBitMemory::GetDWord(uint64 &data)
{
    uint32 c;
    bool r;
    
    r = GetByte(c);
    data = c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    
    return r;
}

// ********************************************************************************************
bool CBitMemory::Get2Bytes(uint32 &data)
{
    uint32 c;
    bool r;
    
    r = GetByte(c);
    data = c;
    r &= GetByte(c);
    data = (data << 8) + c;
    
    return r;
}

// ********************************************************************************************
bool CBitMemory::SkipBytes(uint32 n_bytes)
{
    if(mem_buffer_pos + n_bytes >= mem_buffer_size)
        return false;
    
    mem_buffer_pos += n_bytes;
    
    return true;
}

#endif
