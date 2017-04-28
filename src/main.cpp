/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include <iostream>
#include <string>
#include "params.h"
#include "bit_memory.h"
#include "defs.h"
#include "queues.h"
#include "VCFManager.h"
#include "compression_settings.h"
#include "block_init_compressor.h"
#include "end_compressor.h"
#include "decompressor.h"

using namespace std;

Params params;

int usage();
int usage_compress();
int usage_query();
int compress_input();
int compress_parse_param(int argc, const char *argv[]);
int query(int argc, const char *argv[]);
int query_parse_param(int argc, const char *argv[]);

#ifdef DEVELOPMENT_MODE
int usage_compress_dev();
int usage_query_dev();
int compress_input_dev(int argc, const char *argv[]);
int compress_parse_param_dev(int argc, const char *argv[]);
int query_dev(int argc, const char *argv[]);
int query_parse_param_dev(int argc, const char *argv[]);
#endif

int main(int argc, const char * argv[]) {
    
    if (argc < 2)
        return usage();
    
    if(strcmp(argv[1], "compress") == 0)
    {
        params.task = tcompress;
        if(compress_parse_param(argc, argv) == 1)
            return 1;
        return compress_input();
    }
    else  if(strcmp(argv[1], "view") == 0)
    {
        params.task = tquery;
        return query(argc, argv);
    }
#ifdef DEVELOPMENT_MODE
    else  if(strcmp(argv[1], "compress_dev") == 0)
    {
        
        if(compress_parse_param_dev(argc, argv) == 1)
            return 1;
        return compress_input_dev(argc, argv);
    }
    else  if(strcmp(argv[1], "view_dev") == 0)
    {
        params.task = tquery_dev;
        return query_dev(argc, argv);
    }
#endif
    else
        return usage();
    
}

//--------------------------------------------------------------------------------
// Show execution options
int usage()
{
    cout << "Usage: gtc [option] [arguments] "<< endl;
    cout << "Available options: "<< endl;
    cout << "\tcompress - compress and index VCF/BCF file"<< endl;
    cout << "\tview     - query archive"<< endl;
#ifdef DEVELOPMENT_MODE
    cout << "\tcompress_dev \t- preprocess VCF/BCF file (create BV and IND files) or compress BV+IND files" << endl;
    cout << "\tview_dev \t- query archive and output bit vector with genotypes" << endl;
#endif
    exit (1);
}
int usage_compress()
{
    cout << "Compress VCF/BCF file"<< endl;
    cout << "Usage: gtc compress <options> [file_name]  "<< endl;
    cout << "[file_name]\t\t- input file (a VCF or VCF.GZ file by default)\n"<< endl;
    cout << "Available options (optional): "<< endl;
    cout << "Input: "<< endl;
    cout << "\t-b    \t- input is a BCF file (input is a VCF or VCF.GZ file by default)\t"<< endl;
    cout << "\t-p [x]\t- set ploidy of samples in input VCF to [x] (number >= 1; 2 by default)"<< endl;
    cout << "Output: "<< endl;
    cout << "\t-o [name]\t- set archive name to [name](\"archive\" by default)\t"<< endl;
    cout << "Parameters: "<< endl;
    cout << "\t-t [x]\t- set number of threads to [x] (number >= 1; 2 by default)"<< endl;
    cout << "\t-d [x]\t- set maximum depth to [x] (number >= 0; 0 means no matches; 100 by default)"<< endl;
    cout << "\t-g [x]   \t- [DEV] set number of vector groups [percentage of 1s] to [x] (max: "<< MAX_NUMBER_OF_GROUP <<"; 8 by default)\t"<< endl;
    cout << "\t-hm [x]   \t- [DEV] set n_vec_history for matches to pow(2, [x]) (8 by default, min: " << MATCH_BITS_HUF << ")\t"<< endl;
//    cout << "\t-hc [x]   \t- [DEV] set n_vec_history for copies to pow(2, [x]) (17 by default, min: " << MATCH_BITS_HUF << ")\t"<< endl;
    
    exit (1);
}

int usage_query()
{
    cout << "Output VCF/BCF file"<< endl;
    cout << "Usage: gtc view <options> [archive_name]"<< endl;
    cout << "Available options: "<< endl;
    cout << "Output: "<< endl;
    cout << "\t-o [name]\t- output to a file and set output name to [name] (stdout by default)\t"<< endl;
    cout << "\t-b\t- output a BCF file (output is a VCF file by default)\t"<< endl;
    cout << "\t-C \t- write AC/AN to the INFO field (always set when using -minAC, -maxAC, -minAF or -maxAF)"<< endl;
    cout << "\t-G \t- don't output sample genotypes (only #CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO columns)" <<endl;
    cout << "\t-c [0-9]   set level of compression of the output bcf (number from 0 to 9; 1 by default; 0 means no compression)\t"<< endl;
    cout << "Query: "<< endl;
    cout << "\t-r\t- range in format [chr]:[start]-[end] (for example: -r 14:19000000-19500000). By default all variants are decompressed."<< endl;
    cout << "\t-s\t- sample name(s), separated by comms (for example: -s HG00096,HG00097) OR '@' sign followed by the name of a file with sample name(s) separated by whitespaces (for exaple: -s @file_with_IDs.txt). By default all samples/individuals are decompressed"<< endl;
    cout << "\t-n X \t- process at most X records (by default: all from the selected range)" << endl;
    cout << "Settings: "<< endl;
    cout << "\t-minAC X \t- report only sites with count of alternate alleles among selected samples smaller than or equal to X (default: no limit)" << endl;
    cout << "\t-maxAC X \t- report only sites with count of alternate alleles among selected samples greater than or equal to X" << endl;
    cout << "\t-minAF X \t- report only sites with allele frequency among selected samples greather than or equal to X (X - number between 0 and 1; default: 0)" << endl;
    cout << "\t-maxAF X \t- report only sites with allele frequency among selected samples smaller than or equal to X (X - number between 0 and 1; default: 1)" << endl;
    cout << "\t-m X\t- limit maximum memory usage to remember previous vectors to X MB (no limit by default)\t"<< endl;
    cout << endl;
    exit (1);
}


int compress_input()
{
    VCFManager managerVCF(params);
    if(!managerVCF.CreateIndividualListFile()) return 1;
    uint32_t no_samples = managerVCF.getInputNoSamples();
  
    CompSettings settings(params, no_samples);
    
    CBlockQueue inBlockQueue(max((int) params.n_threads * 2, 8));
    CCompressedBlockQueue compBlockQueue;
    managerVCF.setQueue(&inBlockQueue);
    
    // Distribute blocks to threads, thread initially compresses block and pushes it into compBlockQueue queue
    vector<thread *> workers(params.n_threads, nullptr);
    for(uint32_t i = 0; i < params.n_threads; ++i)
        workers[i] = new thread([&]{
            int id_block = 0;
            unsigned long n_rec;
            unsigned char * ptr  = nullptr, *compressedBlock = nullptr;
            uint32_t * origin_of_copy;
            size_t compressed_size;
            
            vector<int> perm;
            perm.clear();
            perm.resize(no_samples * params.ploidy, 0);
          
            BlockInitCompressor init_compr(&settings);
            while(true)
            {
                
                if(!inBlockQueue.Pop(id_block, ptr, n_rec))
                    break;				// End of blocks data
               
                vector<bool> zeros_only(n_rec, false);
                vector<bool> copies(n_rec, false);
                
                init_compr.SetBlock(n_rec, ptr);
                
                // Permutations
                init_compr.PermuteBlock(perm, true);
                
                // Initial compression
                init_compr.Compress(zeros_only, copies, compressedBlock, compressed_size, origin_of_copy);
                
                compBlockQueue.Push(id_block, compressedBlock, n_rec, compressed_size, perm, zeros_only, copies, origin_of_copy);
        
                delete [] ptr;
            }
        });
    
    if(!managerVCF.ProcessInVCF()) return 1;
    uint64_t no_vec = managerVCF.getNoVec();
    
    managerVCF.CloseFiles();
    
    for(auto p : workers)
    {
        p->join();
        delete p;
    }
    workers.clear();
    std::cout << "All blocks initially compressed." << std::endl;

    // Gathering blocks
    int id_block = 0;
    unsigned long n_rec;
    unsigned char  *compressedBlock = nullptr;
    size_t compressed_size;
    vector<int> perm;
    vector<bool> zeros;
    vector<bool> copies;
    uint32_t * origin_of_copy = nullptr;
    
    EndCompressor endCompressor(&settings, no_vec);
    
    while(compBlockQueue.Pop(id_block, compressedBlock, n_rec, compressed_size, perm, zeros, copies, origin_of_copy))
    {
        endCompressor.AddBlock(id_block, compressedBlock, n_rec, compressed_size, perm, zeros, copies, origin_of_copy);
    
        delete [] compressedBlock;
        delete [] origin_of_copy;
    }
    std::cout << "Final encoding." << std::endl;
    endCompressor.Encode();
    
    endCompressor.storeArchive(params.arch_name.c_str());
    std::cout << "Archive file (" << params.arch_name + ".gtc" << ") created." << std::endl;
    return 0;
}

// Parse the parameters
int compress_parse_param(int argc, const char *argv[])
{
    int i;
    int tmp;
    
    if(argc < 3)
        return usage_compress();
    
    for(i = 2 ; i < argc - 1; ++i)
    {
        if(argv[i][0] != '-')
            break;
        if(strncmp(argv[i], "-b", 2) == 0)
            params.in_type = BCF;
        else if(strncmp(argv[i], "-o", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_compress();
            params.arch_name = string(argv[i]);
        }
        else if(strncmp(argv[i], "-p", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_compress();
            tmp = atoi(argv[i]);
            if(tmp < 1)
                usage_compress();
            
            params.ploidy= tmp;
        }
        else if(strncmp(argv[i], "-d", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_compress();
            tmp = atoi(argv[i]);
            if(tmp < 0)
                usage_compress();
            
            params.max_depth = tmp;
        }
        else if(strncmp(argv[i], "-t", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_compress();
            tmp = atoi(argv[i]);
            if(tmp < 1)
                usage_compress();
            
            params.n_threads = tmp;
        }
        else if(strncmp(argv[i], "-g", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_compress();
            params.ones_ranges = atoi(argv[i]);
            if(params.ones_ranges > MAX_NUMBER_OF_GROUP || params.ones_ranges < 1)
            {
                cout << "wrong value ones_ranges = " << params.ones_ranges << ", (max is " << MAX_NUMBER_OF_GROUP << ", min is 1)" << endl;
                usage_compress();
            }
        }
        else if(strncmp(argv[i], "-hm", 3) == 0)
        {
            i++;
            if(i >= argc)
                return usage_compress();
            if(atoi(argv[i]) < (int) MATCH_BITS_HUF)
            {
                cout << "wrong value hm = " << atoi(argv[i]) << endl;
                usage_compress();
            }
            params.max_bit_size_id_match_pos_diff = atoi(argv[i]); //(1 << atoi(argv[i])) - 1 ;
            
        }
        else if(strncmp(argv[i], "-hc", 3) == 0)
        {
            i++;
            if(i >= argc)
                return usage_compress();
            if(atoi(argv[i]) < (int) MATCH_BITS_HUF)
            {
                cout << "wrong value hc = " << atoi(argv[i]) << endl;
                usage_compress();
            }
            params.max_bit_size_id_copy_pos_diff = atoi(argv[i]); //(1 << atoi(argv[i])) - 1 ;
            
        }
    }
    if(i >= argc)
        return usage_compress();
    
    params.in_file_name = string(argv[i]);
    
    return 0;
}

int query(int argc, const char *argv[])
{
    if(query_parse_param(argc, argv) == 1)
        return 1;
    
    Decompressor decompressor(params); // Load settings and data
    
    decompressor.loadPack();
    decompressor.loadBCF();
    decompressor.setMemoryUsage();
    
    decompressor.initOut();
    decompressor.decompress();

    return 0;
}

// Parse the parameters
int query_parse_param(int argc, const char *argv[])
{
    int i;
    int tmp;
    double tmp_dbl;
    
    if(argc < 3)
        return usage_query();
    
    for(i = 2 ; i < argc - 1; ++i)
    {
        if(argv[i][0] != '-')
        {
            usage_query();
            break;
        }
        if(strncmp(argv[i], "-s", 2) == 0)
        {
            i++;
            params.samples = string(argv[i]);
        }
        else if(strncmp(argv[i], "-n", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query();
            tmp = atoi(argv[i]);
            if(tmp < 0)
                usage_query();
            else
            {
                params.records_to_process = tmp;
            }
        }
        else if(strncmp(argv[i], "-C", 2) == 0)
            params.out_AC_AN = true;
        else if(strncmp(argv[i], "-G", 2) == 0)
            params.out_genotypes = false;
        else if(strncmp(argv[i], "-minAC", 6) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query();
            tmp = atoi(argv[i]);
            if(tmp < 0)
                usage_query();
            params.minAC = tmp;
            params.out_AC_AN = true;
        }
        else if(strncmp(argv[i], "-maxAC", 6) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query();
            tmp = atoi(argv[i]);
            if(tmp < 0)
                usage_query();
            params.maxAC = tmp;
            params.out_AC_AN = true;
        }
        else if(strncmp(argv[i], "-minAF", 6) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query();
            tmp_dbl = atof(argv[i]);
            if(tmp_dbl < 0.0 || tmp_dbl > 1.0)
                usage_query();
            params.minAF = tmp_dbl;
            params.out_AC_AN = true;
        }
        else if(strncmp(argv[i], "-maxAF", 6) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query();
            tmp_dbl = atof(argv[i]);
            if(tmp_dbl < 0.0 || tmp_dbl > 1.0)
                usage_query();
            params.maxAF = tmp_dbl;
            params.out_AC_AN = true;
        }
        else if(strncmp(argv[i], "-r", 2) == 0)
        {
            i++;
            params.range = string(argv[i]);
        }
        else if(strncmp(argv[i], "-b", 2) == 0)
        {
            params.out_type = BCF;
        }
        else if(strncmp(argv[i], "-m", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query();
            tmp = atoi(argv[i]);
            if(tmp < 0)
                usage_query();
            if(!tmp)
                params.MB_memory = false;
            else
                params.max_MB_memory = tmp;
        }
        else if(strncmp(argv[i], "-o", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query();
            params.out_name = string(argv[i]);
        }
        else if(strncmp(argv[i], "-c", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query();
            tmp = atoi(argv[i]);
            if(tmp < 0 || tmp > 9)
                usage_query();
            else
            {
                if(tmp)
                    params.compression_level = argv[i][0];
                else
                    params.compression_level = 'u';
            }
        }
    }
    if(i >= argc)
        return usage_query();
    
    params.arch_name = string(argv[i++]);
    
    return 0;
}

#ifdef DEVELOPMENT_MODE

int compress_input_dev(int argc, const char *argv[])
{
    uint64_t no_vec = 0;
    uint32_t no_samples = 0;
    
   
    if(params.task == tcompress_dev_pre)  // Preprocess only
    {
        VCFManager managerVCF(params);
        if(!managerVCF.CreateIndividualListFile()) return 1;
        no_samples = managerVCF.getInputNoSamples();
        
        CBlockQueue inBlockQueue(max((int) params.n_threads * 2, 8));

        managerVCF.setQueue(&inBlockQueue);
        
        if(!managerVCF.ProcessInVCF()) return 1;
        no_vec = managerVCF.getNoVec();
        
        managerVCF.CloseFiles();
        
        
        uint64_t vec_len = (no_samples *  params.ploidy) / 8 + (((no_samples * params.ploidy) % 8)?1:0);
        uint64_t bv_size = no_vec * vec_len;
        
        char * bv_buf = new char[bv_size];
        uint64_t bv_pos = 0;
        
        int id_block = 0;
        int cnt_block = 0;
        unsigned long n_rec;
        unsigned char * ptr  = nullptr;
        while(inBlockQueue.Pop(id_block, ptr, n_rec))
        {
            assert(id_block == cnt_block++);
            memcpy(bv_buf + bv_pos, ptr, n_rec * vec_len);
            bv_pos += n_rec * vec_len;
            
            delete [] ptr;
        }
        ofstream bv_out(params.arch_name  + ".bv", ios::binary);
        if(!bv_out.is_open())
        {
            cout << "Could not open " << params.arch_name  << ".bv" <<endl;
            exit(1);
        }
        
        if(no_vec/2 > UINT32_MAX)
        {
            cout << "Too many variants - cannot store in 4 bytes.";
            exit(2);
        }
        uint32_t no_var = no_vec/2;
        
        bv_out.write((char *)&no_samples, sizeof(uint32_t));
        bv_out.write((char *)&params.ploidy, sizeof(uint32_t));
        bv_out.write((char *)&no_var, sizeof(uint32_t));
        bv_out.write(bv_buf, bv_pos);
        bv_out.close();
        
        std::cout << "File with bit vectors (" << params.arch_name + ".bv" << ") created." << std::endl;
        
        delete [] bv_buf;
        return 0;
    }
    else //tcompress_dev
    {
        ifstream bv_out(params.arch_name  + ".bv", ios::binary);
        if(!bv_out.is_open())
        {
            cout << "Could not open " << params.arch_name  << ".bv" <<endl;
            exit(1);
        }
        
        bv_out.read((char *)&no_samples, sizeof(uint32_t));
        bv_out.read((char *)&params.ploidy, sizeof(uint32_t));
        bv_out.read((char *)&no_vec, sizeof(uint32_t)); //no_var read
        
        no_vec = no_vec*2;
        uint64_t vec_len = (no_samples *  params.ploidy) / 8 + (((no_samples * params.ploidy) % 8)?1:0);
        
        uint32_t no_vec_in_block = params.var_in_block * 2;
        CompSettings settings(params, no_samples);
        
        CBlockQueue inBlockQueue(max((int) params.n_threads * 2, 8));
        
        CCompressedBlockQueue compBlockQueue;
        
        vector<thread *> workers(params.n_threads, nullptr);
        for(uint32_t i = 0; i < params.n_threads; ++i)
            workers[i] = new thread([&]{
                int id_block = 0;
                unsigned long n_rec;
                unsigned char * ptr  = nullptr, *compressedBlock = nullptr;
                uint32_t * origin_of_copy;
                size_t compressed_size;
                
                vector<int> perm;
                perm.clear();
                perm.resize(no_samples * params.ploidy, 0);
                
                BlockInitCompressor init_compr(&settings);
                while(true)
                {
                    
                    if(!inBlockQueue.Pop(id_block, ptr, n_rec))
                        break;				// End of blocks
                    
                    
                    vector<bool> zeros_only(n_rec, false);
                    vector<bool> copies(n_rec, false);
                    
                    init_compr.SetBlock(n_rec, ptr);
                    
                    // Permutations
                    init_compr.PermuteBlock(perm, true);
                    
                    // Initial comprassion
                    init_compr.Compress(zeros_only, copies, compressedBlock, compressed_size, origin_of_copy);
                    
                    compBlockQueue.Push(id_block, compressedBlock, n_rec, compressed_size, perm, zeros_only, copies, origin_of_copy);
    
                    delete [] ptr;
                }
            });
        
        
        int block_id = 0;
        
        
        uint32_t read_vec = 0;
        
        while(read_vec + no_vec_in_block < no_vec)
        {
            unsigned char * bv_buf = new unsigned char[no_vec_in_block * vec_len]; //deleted by workers
            
            bv_out.read((char *)bv_buf, no_vec_in_block * vec_len);
            
            inBlockQueue.Push(block_id, bv_buf, no_vec_in_block);
            block_id++;
            
            read_vec += no_vec_in_block;
        }
        // Last
        unsigned char * bv_buf = new unsigned char[(no_vec - read_vec) * vec_len]; //deleted by workers
        bv_out.read((char *)bv_buf, (no_vec - read_vec) * vec_len);
        
        inBlockQueue.Push(block_id, bv_buf, (no_vec - read_vec));
        block_id++;
        read_vec += (no_vec - read_vec);
        bv_out.close();
        inBlockQueue.Complete();
        
        for(auto p : workers)
        {
            p->join();
            delete p;
        }
        workers.clear();
        
        
        // Gathering blocks
        int id_block = 0;
        unsigned long n_rec;
        unsigned char  *compressedBlock = nullptr;
        size_t compressed_size;
        vector<int> perm;
        vector<bool> zeros;
        vector<bool> copies;
        uint32_t * origin_of_copy = nullptr;
        
        EndCompressor endCompressor(&settings, no_vec);
        
        while(compBlockQueue.Pop(id_block, compressedBlock, n_rec, compressed_size, perm, zeros, copies, origin_of_copy))
        {
            endCompressor.AddBlock(id_block, compressedBlock, n_rec, compressed_size, perm, zeros, copies, origin_of_copy);
            
            delete [] compressedBlock;
            delete [] origin_of_copy;
        }
        
        endCompressor.Encode();
        endCompressor.storeArchive(params.arch_name.c_str());
    }
    
    return 0;
}


// Parse the parameters query_dev
int query_parse_param_dev(int argc, const char *argv[])
{
    int i;
    int tmp;
    
    if(argc < 3)
        return usage_query_dev();
    
    params.out_type = BV;
    for(i = 2 ; i < argc - 1; ++i)
    {
        if(argv[i][0] != '-')
            break;
        else if(strncmp(argv[i], "-m", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query_dev();
            tmp = atoi(argv[i]);
            if(tmp < 0)
                usage_query_dev();
            if(!tmp)
                params.MB_memory = false;
            else
                params.max_MB_memory = tmp;
        }
        else if(strncmp(argv[i], "-o", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query_dev();
            params.out_name = string(argv[i]);
        }
        else if(strncmp(argv[i], "-t", 2) == 0)
        {
            params.out_type = TXT_BV;
        }
        else if(strncmp(argv[i], "-i", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query_dev();
            tmp = atoi(argv[i]);
            if(tmp < 0)
                usage_query_dev();
            params.sample_to_dec = tmp;
            params.dec_single_sample = true;
        }
        else if(strncmp(argv[i], "-j", 2) == 0)
        {
            i++;
            if(i >= argc)
                return usage_query_dev();
            tmp = atoi(argv[i]);
            if(tmp < 0)
                usage_query_dev();
            params.var_to_dec = tmp;
            params.dec_single_var = true;
        }
    }
    if(i >= argc)
        return usage_query_dev();
    
    params.arch_name = string(argv[i++]);
    return 0;
}

int query_dev(int argc, const char *argv[])
{
    if(query_parse_param_dev(argc, argv) == 1)
        return 1;
    
    Decompressor decompressor(params); // Load settings and data
    
    decompressor.loadPack();
    decompressor.setMemoryUsage();
    
    if(params.dec_single_var == false && params.dec_single_sample == false)
    {
        decompressor.initBVOut(true);
        decompressor.getBV();
    }
    else if(params.dec_single_var == true && params.dec_single_sample == true)
    {
        decompressor.initBVOut(false);
        decompressor.getVarSampleBV(params.var_to_dec, params.sample_to_dec);
    }
    else if(params.dec_single_var == true && params.dec_single_sample == false)
    {
        decompressor.initBVOut(false);
        decompressor.getVarBV(params.var_to_dec);
    }
    else //if(params.dec_single_var == false && params.dec_single_sample == true)
    {
        decompressor.initBVOut(false);
        decompressor.getSampleBV(params.sample_to_dec);
    }
  
    return 0;
}

// Parse the parameters
int compress_parse_param_dev(int argc, const char *argv[])
{
    int i;
    int tmp;
    
    if(argc < 4)
        return usage_compress_dev();
    
    if(argv[2][0] == 'p')
    {
        params.task = tcompress_dev_pre;
        params.mode = 'p';
        params.preprocessVCFonly = true;
        for(i = 3 ; i < argc - 1; ++i)
        {
            if(argv[i][0] != '-')
                break;
            if(strncmp(argv[i], "-b", 2) == 0)
                params.in_type = BCF;
            else if(strncmp(argv[i], "-o", 2) == 0)
            {
                i++;
                if(i >= argc)
                    return usage_compress_dev();
                params.arch_name = string(argv[i]);
            }
            else if(strncmp(argv[i], "-p", 2) == 0)
            {
                i++;
                if(i >= argc)
                    return usage_compress_dev();
                tmp = atoi(argv[i]);
                if(tmp < 1)
                    usage_compress_dev();
                
                params.ploidy= tmp;
            }
        }
        if(i >= argc)
            return usage_compress_dev();
        
        params.in_file_name = string(argv[i]);
        
    }
    else if(argv[2][0] == 'c')
    {
        params.task = tcompress_dev;
        
        params.mode = 'c';
        params.input_is_bit_vector = true;
        for(i = 3 ; i < argc - 2; ++i)
        {
            if(argv[i][0] != '-')
                break;
            
            if(strncmp(argv[i], "-d", 2) == 0)
            {
                 i++;
                 if(i >= argc)
                     return usage_compress_dev();
                 tmp = atoi(argv[i]);
                 if(tmp < 0)
                     usage_compress_dev();
                 
                 params.max_depth = tmp;
             }
             else if(strncmp(argv[i], "-t", 2) == 0)
             {
                 i++;
                 if(i >= argc)
                     return usage_compress();
                 tmp = atoi(argv[i]);
                 if(tmp < 1)
                     usage_compress();
                 
                 params.n_threads = tmp;
             }
             else if(strncmp(argv[i], "-g", 2) == 0)
             {
                 i++;
                 if(i >= argc)
                     return usage_compress_dev();
                 params.ones_ranges = atoi(argv[i]);
                 if(params.ones_ranges > MAX_NUMBER_OF_GROUP || params.ones_ranges < 1)
                 {
                     cout << "wrong value ones_ranges = " << params.ones_ranges << ", (max is " << MAX_NUMBER_OF_GROUP << ", min is 1)" << endl;
                     usage_compress_dev();
                 }
             }
             else if(strncmp(argv[i], "-hm", 3) == 0)
             {
                 i++;
                 if(i >= argc)
                     return usage_compress_dev();
                 if(atoi(argv[i]) < MATCH_BITS_HUF)
                 {
                     cout << "wrong value hm = " << atoi(argv[i]) << endl;
                     usage_compress_dev();
                 }
                 params.max_bit_size_id_match_pos_diff = atoi(argv[i]); //(1 << atoi(argv[i])) - 1 ;
                 
             }
             else if(strncmp(argv[i], "-hc", 3) == 0)
             {
                 i++;
                 if(i >= argc)
                     return usage_compress_dev();
                 if(atoi(argv[i]) < MATCH_BITS_HUF)
                 {
                     cout << "wrong value hc = " << atoi(argv[i]) << endl;
                     usage_compress_dev();
                 }
                 params.max_bit_size_id_copy_pos_diff = atoi(argv[i]); //(1 << atoi(argv[i])) - 1 ;
                 
             }
        }
        if(i >= argc)
            return usage_compress_dev();
        
        params.arch_name = string(argv[i]);
    }
    else
        return usage_compress_dev();
    
    return 0;
}

int usage_query_dev()
{
    cout << "[DEV, ALTERNATIVE USAGE] \nQuery archive (only *.gtc file is read/needed) \nand output bit vector with genotypes (no information about position, alleles or sample ID). \nIt is possible to query for all collection (default), i-th sample and/or j-th variant." << endl<< endl;
    cout << "Usage: gtc view_dev <options> [archive_name]"<< endl;
    cout << "[archive_name] - name of the archive (archive_name.gtc file is needed)"<< endl;
    cout << endl;
    cout << "Available options: "<< endl;
    
    cout << "Output: "<< endl;
    cout << "\t-o [name]\t- output to a file and set output name to [name] (stdout by default)\t"<< endl;
    cout << "\t-t [name]\t- output text bit/byte vector to standard output (not by default)\t"<< endl;
    
    cout << "View/query: "<< endl;
    cout << "\t-i I output bit vector with genotype of I-th sample (all by default; 0-based); if used with -j one byte for each haplotype is outputted (instead of bit)"<< endl;
    cout << "\t-j J output bit vector with genotypes at J-th variant site (all by default; 0-based); if used with -i one byte for each haplotype is outputted (instead of bit)"<< endl;
    
    cout << "Settings: "<< endl;
    cout << "\t-m X\t- limit maximum memory usage to remember previous vectors to X MB (no limit by default)\t"<< endl;
    exit (1);
}

int usage_compress_dev()
{
    cout << "[DEV, ALTERNATIVE USAGE] \nOption p to preprocess VCF/BCF file to create bit vector with genotypes (*.bv; no information about position, alleles or sample ID; header 12 bytes: no_samples, ploidy, no_var) and matching file with all samples IDs (*.ind). It also creates a *.bcf file with variants description.\nOption c to compress [name].bv file, creates [name].gtc " << endl<< endl;
    
    cout << "Usage: gtc compress_dev [mode] <options> [in_name]"<< endl;
    cout << "[mode]\t- mode: p for preprocessing of VCF/BCF [in_name] file;\n      \t        c for compression of intermediate bv file ([in_name].bv)"<< endl;
    cout << "[in_name]\t- input file: single VCF file for mode p; name of preprocessed bit vector file ([in_name].bv) for mode c\n"<< endl<< endl;
    cout << "Mode p options (optional): "<< endl;
    cout << "\t-b    \t- input is a BCF file (input is a VCF or VCF.GZ file by default)\t"<< endl;
    cout << "\t-p [x]\t- set ploidy of samples in input VCF/BCF to [x] (number >= 1; 2 by default)"<< endl;
    cout << "\t-o [name]\t- set output files names to [name] (\"archive\" by default)\t"<< endl<< endl;
    
    cout << "Mode c options (optional): "<< endl;
    cout << "\t-d [x]\t- set maximum depth to [x] (number >= 0; 0 means no matches; 100 by default)"<< endl;
    cout << "\t-g [x]   \t- [DEV] set number of vector groups [percentage of 1s] to [x] (max: "<< MAX_NUMBER_OF_GROUP <<"; 32 by default)\t"<< endl;
    cout << "\t-hm [x]   \t- [DEV] set n_vec_history for matches to pow(2, [x]) (9 by default, min: " << MATCH_BITS_HUF << ")\t"<< endl;
    cout << "\t-hc [x]   \t- [DEV] set n_vec_history for copies to pow(2, [x]) (17 by default, min: " << MATCH_BITS_HUF << ")\t"<< endl;
    cout << "\t-t [x]\t- set number of threads to [x] (number >= 1; 8 by default)"<< endl;
    
    exit (1);
}

#endif

