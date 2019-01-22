/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include "defs.h"
#include "samples.h"
#include "htslib/vcf.h"
#include <iostream>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <algorithm>

int Samples::loadSamples(std::string filename)
{
    std::string ind_name;
    std::ifstream fi(filename);
    
    if(!fi.is_open())
        return 1;
    
    uint32_t c = 0;
    while(fi >> ind_name)
    {
        if(!whichIndMap.insert ( std::pair<std::string,uint32_t>(ind_name, c++) ).second)
        {
            std::cout << "Error! Two individuals with the same name!\n";
            return 2;
        }
        
    }
    fi.close();
    
    return 0;
}


uint32_t Samples::getWhich(std::string nm)
{
    auto it = whichIndMap.find(nm);
    if(it != whichIndMap.end())
             return it->second;
    else
    {
        std::cout << "There is no sample " << nm << " in the set\n";
        exit(1);
        
    }

}

int Samples::setAllSamples(bcf_hdr_t * hdr, std::string filename, bool out_genotypes)
{
    std::string ind_name;
    std::ifstream fi(filename);
    
    if(!fi.is_open())
        return 1;
    
    uint32_t c = 0;
    while(fi >> ind_name)
    {
        if(!whichIndMap.insert ( std::pair<std::string,uint32_t>(ind_name, c++) ).second)
        {
            std::cout << "Error! Two individuals with the same name!\n";
            return 2;
        }
        if(out_genotypes)
            bcf_hdr_add_sample(hdr, ind_name.c_str());
        
    }
    fi.close();
    
    no_samples = (uint32_t)whichIndMap.size();
    return 0;

}

uint32_t * Samples::setSamples(bcf_hdr_t * hdr, const std::string & samples, bool out_genotypes)
{
    uint32_t * smplIDs = nullptr;
    long size = 0;
    
    no_samples = 0;
    if(samples[0] == '@')
    {
        std::ifstream in_samples(samples.substr(1));
        if(!in_samples.is_open())
        {
            std::cout << "Error. Cannot open " << samples.substr(1)<< " file with samples.\n";
            exit(1);
        }
        
        uint32_t i = 0;
        std::string item;
        
        while (in_samples >> item) {
            size++;
        }
        in_samples.clear();
        in_samples.seekg(0);
        
        smplIDs = new uint32_t[size];
        while (in_samples >> item) {
            if(out_genotypes)
                bcf_hdr_add_sample(hdr, item.c_str());
            smplIDs[i++] = getWhich(item);
            no_samples++;
        }

    }
    else
    {
        size = std::count(samples.begin(), samples.end(), ',') + 1;
        
        smplIDs = new uint32_t[size];
        char delim = ',';
        uint32_t i = 0;
        std::stringstream ss(samples);
        std::string item;
        while (getline(ss, item, delim)) {
            if(out_genotypes)
                bcf_hdr_add_sample(hdr, item.c_str());
            smplIDs[i++] = getWhich(item);
            no_samples++;
        }

    }
    
    return smplIDs;
}
