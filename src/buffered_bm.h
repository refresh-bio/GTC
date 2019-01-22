/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */
#ifndef buffered_bm_h
#define buffered_bm_h

#include "defs.h"
#include <iostream>
#include "bit_memory.h"
#include "huffman.h"

class CBufferedBitMemory
{
    CBitMemory * bm = NULL;
    uint32_t buffer = 0;
    
    uint8_t no_bits;
    uint32_t temp;
    
public:
    uint16_t left = 0;
    CBufferedBitMemory();
    void setBitMemory(CBitMemory * _bm);
    uint32_t getBits(uint32_t no_bits);
    void getBitsAndDiscard(uint32_t no_bits);
    
    bool SetPos(int64 pos);
    
    int32_t decodeFastLut(CHuffman * huff); //first get max len
    int32_t decodeFastLutStat(CHuffman * huff, uint32_t & no_read_bits);
    int32_t decodeFast(CHuffman * huff); //first get min len
    int32_t decodeFastStat(CHuffman * huff, uint32_t &no_read_bits);
    
    void setBuffer(uint32_t _buffer, uint16_t _left);
    void getBuffer(uint32_t & _buffer, uint16_t & _left);
};

#endif /* buffered_bm_h */
