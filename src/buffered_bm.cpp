/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include "buffered_bm.h"

CBufferedBitMemory::CBufferedBitMemory()
{
    bm = NULL;
    left = 0;
    buffer = 0;
};

void CBufferedBitMemory::setBitMemory(CBitMemory * _bm)
{
    bm = _bm;
}

bool CBufferedBitMemory::SetPos(int64 pos)
{
    left = 0;
    buffer = 0;
    return bm->SetPos(pos);
}

int32_t CBufferedBitMemory::decodeFastLut(CHuffman * huff)
{
    int32_t result;
    uint32_t toGet;
    uint32_t max_lut_len = huff->max_lut_len;
    
    if(left)
    {
        if(left < max_lut_len)
        {
            buffer = buffer & bm->n_bit_mask[left];
            toGet = max_lut_len - left;
            temp = buffer << toGet;
            bm->GetBits(buffer, toGet);
            buffer |= temp;
            
            result = huff->DecodeFastLut(buffer, no_bits);
            while (result < 0)
            {
                bm->GetBit(buffer);
                result = huff->Decode(buffer);
            }
            left = max_lut_len - no_bits;

            return result;
        }
        else
        {
            buffer = buffer & bm->n_bit_mask[left];
            left -= max_lut_len;
            temp = buffer >> left;
            
            result = huff->DecodeFastLut(temp, no_bits);
            
            while (result < 0)
            {
                if(left)
                {
                    temp = (buffer >> (left--  - 1) ) & 1;
                    result = huff->Decode(temp);
                }
                else
                {
                    bm->GetBit(buffer);
                    result = huff->Decode(buffer);
                }
            }
            left = max_lut_len - no_bits + left;
            return result;
            
        }
    }
    else
    {
        bm->GetBits(buffer, max_lut_len);
        result = huff->DecodeFastLut(buffer, no_bits);
        while (result < 0)
        {
            bm->GetBit(buffer);
            result = huff->Decode(buffer);
        }
        left = max_lut_len - no_bits;
        return result;
    }
}

int32_t CBufferedBitMemory::decodeFastLutStat(CHuffman * huff, uint32_t & no_read_bits)
{
    no_read_bits = 0;
    int32_t result;
    uint32_t toGet;
    uint32_t max_lut_len = huff->max_lut_len;
    
    if(left)
    {
        if(left < max_lut_len)
        {
            buffer = buffer & bm->n_bit_mask[left];
            toGet = max_lut_len - left;
            temp = buffer << toGet;
            bm->GetBits(buffer, toGet);
            buffer |= temp;
            
            result = huff->DecodeFastLut(buffer, no_bits);
            no_read_bits += no_bits;
            while (result < 0)
            {
                bm->GetBit(buffer);
                result = huff->Decode(buffer);
                no_read_bits++;
            }
            left = max_lut_len - no_bits;
            
            return result;
        }
        else
        {
            buffer = buffer & bm->n_bit_mask[left];
            left -= max_lut_len;
            temp = buffer >> left;
            
            result = huff->DecodeFastLut(temp, no_bits);
            no_read_bits += no_bits;
            
            while (result < 0)
            {
                if(left)
                {
                    temp = (buffer >> (left--  - 1) ) & 1;
                    result = huff->Decode(temp);
                    no_read_bits++;
                }
                else
                {
                    bm->GetBit(buffer);
                    result = huff->Decode(buffer);
                    no_read_bits++;
                }
            }
            left = max_lut_len - no_bits + left;
            return result;
        }
    }
    else
    {
        bm->GetBits(buffer, max_lut_len);
        result = huff->DecodeFastLut(buffer, no_bits);
        no_read_bits += no_bits;
        while (result < 0)
        {
            bm->GetBit(buffer);
            result = huff->Decode(buffer);
        }
        left = max_lut_len - no_bits;
        return result;
    }
}

int32_t CBufferedBitMemory::decodeFast(CHuffman * huff)
{
    uint32_t min_len = huff->min_len;
    
    if(!huff->max_len)
        return 0;
    int32_t result;
    
    if(left)
    {
        if(left > min_len)
        {
            buffer = buffer & bm->n_bit_mask[left];
            left = left - min_len;
            temp = buffer >> (left);
                
            result = huff->DecodeFast(temp);
            while (result < 0)
            {
                if(left)
                {
                    temp = (buffer >> (left--  - 1) ) & 1;
                    result = huff->Decode(temp);
                }
                else
                {
                    bm->GetBit(buffer);
                    result = huff->Decode(buffer);
                }
            }
                
            return result;
        }
        else
        {
            buffer = buffer & bm->n_bit_mask[left];
            temp = buffer << (min_len - left);
            bm->GetBits(buffer, min_len - left);
            left = 0;
            buffer = buffer | temp;

            result = huff->DecodeFast(buffer);
            while (result < 0)
            {
                bm->GetBit(buffer);
                result = huff->Decode(buffer);
            }
                
            return result;
        }
    }
    else
    {
        bm->GetBits(buffer, min_len);
        result = huff->DecodeFast(buffer);
        while (result < 0)
        {
            bm->GetBit(buffer);
            result = huff->Decode(buffer);
        }

        return result;
    }
}

int32_t CBufferedBitMemory::decodeFastStat(CHuffman * huff, uint32_t &no_read_bits)
{
    no_read_bits = 0;
    uint32_t min_len = huff->min_len;
    
    if(!huff->max_len)
        return 0;
    int32_t result;
    
    no_read_bits = min_len;
    
    if(left)
    {
        if(left > min_len)
        {
            buffer = buffer & bm->n_bit_mask[left];
            left = left - min_len;
            temp = buffer >> (left);
            
            result = huff->DecodeFast(temp);
            while (result < 0)
            {
                if(left)
                {
                    temp = (buffer >> (left--  - 1) ) & 1;
                    result = huff->Decode(temp);
                    no_read_bits++;
                }
                else
                {
                    bm->GetBit(buffer);
                    result = huff->Decode(buffer);
                    no_read_bits++;
                }
            }
            
            return result;
        }
        else
        {
            buffer = buffer & bm->n_bit_mask[left];
            temp = buffer << (min_len - left);
            bm->GetBits(buffer, min_len - left);
            left = 0;
            buffer = buffer | temp;
            
            result = huff->DecodeFast(buffer);
            while (result < 0)
            {
                bm->GetBit(buffer);
                result = huff->Decode(buffer);
                no_read_bits++;
            }
            
            return result;
            
        }
    }
    else
    {
        bm->GetBits(buffer, min_len);
        result = huff->DecodeFast(buffer);
        while (result < 0)
        {
            bm->GetBit(buffer);
            result = huff->Decode(buffer);
            no_read_bits++;
        }

        return result;
    }
}

uint32_t CBufferedBitMemory::getBits(uint32_t no_bits)
{
    uint32_t result = 0, toGet;
    
    if( left > no_bits)
    {
        buffer = buffer & bm->n_bit_mask[left];
        result = buffer >> (left - no_bits);
        left = left - no_bits;

        return result;
    }
    else
    {
        if(left)
        {
            buffer = buffer & bm->n_bit_mask[left];
            toGet = no_bits - left;
            result = buffer << toGet;
            bm->GetBits(temp, toGet);
            result = result | temp;
            left = 0;

            return result;
        }
        else
        {
            bm->GetBits(result, no_bits);

            return result;
        }
    }
}

void CBufferedBitMemory::getBitsAndDiscard(uint32_t no_bits)
{
    if(no_bits > left)
    {
        bm->GetBitsAndDiscard(no_bits - left);
        left = 0;
    }
    else
        left -= no_bits;
}

void CBufferedBitMemory::getBuffer(uint32_t & _buffer, uint16_t & _left)
{
    _buffer = buffer;
    _left = left;
}

void CBufferedBitMemory::setBuffer(uint32_t _buffer, uint16_t _left)
{
    buffer = _buffer;
    left = _left;
}
