/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#ifndef _QUEUES_H
#define _QUEUES_H

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <map>
#include <list>
#include <stack>
#include <tuple>
#include <set>

using namespace std;

// ********************************************************************************
class CBlockQueue
{
    typedef struct block_tag
    {
        int block_id;
        unsigned char *ptr;
        size_t n_recs;
        
        
        block_tag(int _block_id, unsigned char* _ptr, size_t _n_recs) :
        block_id(_block_id), ptr(_ptr), n_recs(_n_recs)
        {}
    } block_t;
    
    queue<block_t> q_blocks;
    
    bool eoq_flag;
    int capacity;
    
    mutex mtx;
    condition_variable cv_pop, cv_push;
    
public:
    CBlockQueue() : eoq_flag(false), capacity(0)
    {}
    
    CBlockQueue(int _capacity) : eoq_flag(false), capacity(_capacity)
    {}
    
    ~CBlockQueue()
    {}
    
    void Push(int id_block, unsigned char *ptr, size_t n_recs)
    {
        unique_lock<std::mutex> lck(mtx);
        cv_push.wait(lck, [this] {return (int) q_blocks.size() < capacity;});
        
        q_blocks.push(block_t(id_block, ptr, n_recs));
        
        cv_pop.notify_all();
    }
    
    bool Pop(int &id_block, unsigned char *&ptr, size_t &n_recs)
    {
        unique_lock<std::mutex> lck(mtx);
        cv_pop.wait(lck, [this] {return !q_blocks.empty() || eoq_flag; });
        
        if (eoq_flag && q_blocks.empty())
            return false;
        
        id_block = q_blocks.front().block_id;
        ptr		 = q_blocks.front().ptr;
        n_recs   = q_blocks.front().n_recs;
        
        q_blocks.pop();
        
        cv_push.notify_all();
        
        return true;
    }
    
    void Complete()
    {
        unique_lock<std::mutex> lck(mtx);
        
        eoq_flag = true;
        
        cv_pop.notify_all();
    }
};

// ********************************************************************************
class CCompressedBlockQueue
{
    typedef struct compressed_block_tag
    {
        int block_id;
        unsigned char *ptr;
        size_t n_recs;
        size_t compressed_size;
        vector<int> perm;
        vector<bool> zeros;
        vector<bool> copies;
        uint32_t * origin_of_copy;
        
        compressed_block_tag(
                             int _block_id, unsigned char *_ptr, size_t _n_recs, size_t _compressed_size, vector<int> &_perm, vector<bool> &_zeros, vector<bool> &_copies, uint32_t * _origin_of_copy) :
        block_id(_block_id), ptr(_ptr), n_recs(_n_recs), compressed_size(_compressed_size), perm(_perm), zeros(_zeros), copies(_copies), origin_of_copy(_origin_of_copy)
        {}
        
        friend bool operator<(const compressed_block_tag &x, const compressed_block_tag &y)
        {
            return x.block_id < y.block_id;
        }
    } compressed_block_t;
    
    set<compressed_block_t> s_blocks;
    
    mutex mtx;
    
public:
    CCompressedBlockQueue()
    {}
    
    ~CCompressedBlockQueue()
    {}
    
    void Push(int id_block, unsigned char *ptr, size_t n_recs, size_t compressed_size, vector<int> &perm, vector<bool> &zeros, vector<bool> &copies, uint32_t * _origin_of_copy)
    {
        lock_guard<std::mutex> lck(mtx);
        
        s_blocks.insert(compressed_block_t(id_block, ptr, n_recs, compressed_size, perm, zeros, copies, _origin_of_copy));
    }
    
    bool Pop(int &id_block, unsigned char *&ptr, size_t &n_recs, size_t &compressed_size, vector<int> &perm, vector<bool> &zeros, vector<bool> &copies, uint32_t *& origin_of_copy)
    {
        unique_lock<std::mutex> lck(mtx);
        
        if (s_blocks.empty())
            return false;
        
        auto x = s_blocks.begin();
        
        id_block = x->block_id;
        ptr = x->ptr;
        n_recs = x->n_recs;
        compressed_size = x->compressed_size;
        perm = x->perm;
        zeros = x->zeros;
        copies = x->copies;
        origin_of_copy = x->origin_of_copy;
        
        s_blocks.erase(s_blocks.begin());
        
        return true;
    }
};

#endif
