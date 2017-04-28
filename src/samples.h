/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#ifndef samples_h
#define samples_h

#include <iostream>
#include <map>
#include <cstring>
#include "htslib/vcf.h"

class Samples{
     std::map<std::string, uint32_t> whichIndMap;
public:
    uint32_t no_samples;
    
    int loadSamples(std::string filename);
    uint32_t getWhich(std::string nm);
    int setAllSamples(bcf_hdr_t * hdr, std::string filename, bool out_genotypes);
    uint32_t * setSamples(bcf_hdr_t * hdr, const std::string  & samples, bool out_genotypes);
};


#endif /* samples_h */
