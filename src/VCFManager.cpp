/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include "VCFManager.h"

bool VCFManager::OpenInVCF()
{
    if(in_open)
    {
        std::cout << "Input VCF/BCF already opened." << std::endl;
        return false;
    }
    if(in_type == VCF)
    {
        in_vcf = hts_open(in_vcf_name.c_str(), "r");
    }
    else //if(in_type == BCF)
    {
        in_vcf = hts_open(in_vcf_name.c_str(), "rb");
    }
    
    if(!in_vcf)
    {
        std::cout << "could not open " << in_vcf_name << " file" << std::endl;
        return false;
    }
    
    in_open = true;
    return true;
}

bool VCFManager::OpenOutVCF()
{
    if(out_open)
    {
        std::cout << "Output VCF/BCF already opened." << std::endl;
        return false;
    }
    
    if(out_type == VCF)
    {
        out_vcf = hts_open(out_vcf_name.c_str(),"w");
    }
    else
    {
        out_vcf = hts_open(out_vcf_name.c_str(),"wb");
    }
    
    if(!out_vcf)
    {
        std::cout << "Could not open " << out_vcf_name << " file" << std::endl;
        return false;
    }
    out_open = true;
    return true;
}

// Get number of samples in VCF
void VCFManager::ReadInHdr()
{
    in_hdr = bcf_hdr_read(in_vcf);
    no_samples = bcf_hdr_nsamples(in_hdr);
    in_hdr_read = true;
}

uint32_t VCFManager::getInputNoSamples() {
    
    if(!in_open)
        if(!OpenInVCF())
            return 0;
    if(!in_hdr_read)
        ReadInHdr();
    
    return no_samples;
}

// Create file with list of sample names
bool VCFManager::CreateIndividualListFile()
{
    if(!in_open)
        if(!OpenInVCF())
            return false;
    if(!in_hdr_read)
        ReadInHdr();
    
    std::ofstream ind_file(arch_name + ".ind");
    if(ind_file)
    {
        for (uint32_t i=0; i< no_samples; i++)
        {
            ind_file << in_hdr->samples[i] << std::endl;
        }
        std::cout << "File with list of samples (" << arch_name + ".ind" << ") created." << std::endl;
        
        ind_file.close();
        return true;
    }
    else
    {
        std::cout << "Could not open "<< arch_name + ".ind" << "file with list of samples."  << std::endl;
        return false;
    }
    return true;
}

void VCFManager::setBitVector()
{
    if(no_samples == 0)
        no_samples = getInputNoSamples();
    
    vec_len =  (no_samples *  ploidy) / 8 + (((no_samples * ploidy) % 8)?1:0) ;
    
    block_max_size = vec_len * no_vec_in_block + 1;
    bv.Create(block_max_size);
    block_id = 0;
    vec_read_in_block = 0;
}

// Splits multiple alleles sites, reads genotypes, creates blocks of bytes to process, fills out [archive_name].bcf file
bool VCFManager::ProcessInVCF()
{
    if(!in_open)
        if(!OpenInVCF())
            return false;
    if(!in_hdr_read)
        ReadInHdr();

    
    if(!out_open)
        if(!OpenOutVCF())
            return false;
    setBitVector();

    no_vec = 0;
    // Remove info fields
    bcf_hdr_remove(in_hdr,BCF_HL_INFO, NULL);
    
    // Set hdr_out without samples/individuals
    out_hdr = bcf_hdr_subset(in_hdr, 0, NULL, NULL);
    bcf_hdr_sync(out_hdr);
    
    // Add description ot "alt in other line" ALT to out_hdr
    bcf_hdr_append(out_hdr, "##INFO=<ID=_row,Number=1,Type=Integer,Description=\"row number\">");
    bcf_hdr_sync(out_hdr);
    
    // Write header to file
    bcf_hdr_write(out_vcf, out_hdr);
    
    int32_t tmpi = 0;
    char * alt_desc = new char[256];
    bcf1_t *rec    = bcf_init1();
    bcf1_t *new_rec = bcf_init1();
    
    while ( bcf_read1(in_vcf, in_hdr, rec)>=0 )
    {
        if(rec->errcode)
        {
            std::cout << "Repair VCF file\n";
            exit(9);
        }
        bcf_unpack(rec, BCF_UN_ALL);  // Unpack all in record
        if(rec->d.fmt->n !=  (int) ploidy)
        {
            std::cout << "Wrong ploidy (not equal to " << ploidy << ") for record at position " << rec->pos <<".\n";
            std::cout << "Repair VCF file OR set correct ploidy using -p option\n";
            exit(9);
        }
        
        // Set CHROM field
        new_rec->rid = rec->rid;
        // Set POS field
        new_rec->pos = rec->pos;
        // Set QUAL field
        new_rec->qual = 0;
        
        if( tmpi%100000 == 0)
            std::cout << tmpi << " variants preprocessed\n";

        if(rec->n_allele > 2)
        {
            //if "ALT,<M>", do not change line(as VCF was already altered)
            if(rec->n_allele==3 && strcmp(rec->d.allele[2], "<M>") == 0)
            {
                // Check if alt_desc size is enough for alleles
                uint64_t allele_size = strlen(rec->d.allele[0])+strlen(rec->d.allele[1])+5;
                if(allele_size > 256)
                {
                    delete [] alt_desc;
                    alt_desc = new char[allele_size];
                }
                
                strcpy(alt_desc, rec->d.allele[0]);
                strcat(alt_desc, ",");
                strcat(alt_desc, rec->d.allele[1]);
                strcat(alt_desc, ",");
                strcat(alt_desc, rec->d.allele[2]);
                
                bcf_update_alleles_str(out_hdr, new_rec, alt_desc);
                
                addGTtoBitVector(in_hdr, rec);
                
                bcf_update_info_int32(out_hdr, new_rec, "_row", &tmpi, 1);
                tmpi++;
                
                bcf_write1(out_vcf, out_hdr, new_rec);
                bcf_clear(new_rec);
                bcf_clear(rec);
            }
            else  // Break multi alleles into several lines/sites
            {
                
                for(int a=1; a < rec->n_allele; a++)  //create one line for each single allele
                {
                                      
                    // Check if alt_desc size is enough for alleles
                    uint64_t allele_size = strlen(rec->d.allele[0])+strlen(rec->d.allele[a])+5;
                    if(allele_size > 256)
                    {
                        delete [] alt_desc;
                        alt_desc = new char[allele_size];
                    }
                    
                    // Cut unnecessary letters at the end
                    size_t len_ref = strlen(rec->d.allele[0]), len_alt = strlen(rec->d.allele[a]);
                    if(len_ref > 1  && len_alt > 1)
                    {
                        while((rec->d.allele[0][len_ref-1] == rec->d.allele[a][len_alt-1] ) && len_ref > 1 && len_alt > 1)
                        {len_ref--; len_alt--;}
                    }
                    
                    strcpy(alt_desc, "\0");
                    strncat(alt_desc, rec->d.allele[0], len_ref);
                    strcat(alt_desc, ",");
                    strncat(alt_desc, rec->d.allele[a], len_alt);
                    strcat(alt_desc, ",<M>");
                    
                    // New record
                    bcf_update_alleles_str(out_hdr, new_rec, alt_desc);
                    
                    int *gt_arr = NULL, ngt_arr = 0;
                    bcf_get_genotypes(in_hdr, rec, &gt_arr, &ngt_arr);
                    for(int i = 0; i < ngt_arr; i++)
                    {   // gt_arr needed to create bit vectors
                        if(bcf_gt_allele(gt_arr[i]) != 0)
                        {
                            if(bcf_gt_allele(gt_arr[i]) == a)
                                gt_arr[i] = bcf_gt_phased(1);
                            else if(bcf_gt_is_missing(gt_arr[i]))
                                ;
                            else
                                gt_arr[i] = bcf_gt_phased(2);
                        }
                    }
                    
                    // New record
                    bcf_update_genotypes(in_hdr, new_rec, gt_arr, bcf_hdr_nsamples(in_hdr)*ploidy);
                    
                    
                    // Add genotypes from site annotation to two vectors in bit memory vec
                    addGTtoBitVector(in_hdr, new_rec);
                    
                    bcf_update_info_int32(out_hdr, new_rec, "_row", &tmpi, 1) ;
                    
                    // New record
                    bcf_subset(out_hdr, new_rec, 0, 0); // Do not write genotypes
                    
                    tmpi++;
                    
                    // New record
                    bcf_write1(out_vcf, out_hdr, new_rec);
                    
                    free(gt_arr);
                    bcf_clear(new_rec);
                    
                    // Set CHROM field
                    new_rec->rid = rec->rid;
                    // Set POS field
                    new_rec->pos = rec->pos;
                    // Set QUAL field
                    new_rec->qual = 0;
                    
                }
                bcf_clear(rec);
            }
        }
        else // Do not change line (single allele)
        {
            // Check if alt_desc size is enough for alleles
            uint64_t allele_size = strlen(rec->d.allele[0])+strlen(rec->d.allele[1])+5;
            if(allele_size > 256)
            {
                delete [] alt_desc;
                alt_desc = new char[allele_size];
            }
            
            strcpy(alt_desc, rec->d.allele[0]);
            strcat(alt_desc, ",");
            strcat(alt_desc, rec->d.allele[1]);
            
            bcf_update_alleles_str(out_hdr, new_rec, alt_desc);
            
            addGTtoBitVector(in_hdr, rec);
            
            bcf_update_info_int32(out_hdr, new_rec, "_row", &tmpi, 1);
            tmpi++;
            
            bcf_write1(out_vcf, out_hdr, new_rec);
            bcf_clear(new_rec);
            bcf_clear(rec);
        }
        
    }
    
    std::cout << "All variants preprocessed\n";

    // Last pack (may be smaller than block size
    if(vec_read_in_block)
    {
       bv.TakeOwnership();
       queue->Push(block_id, bv.mem_buffer, vec_read_in_block);
      
       no_vec = no_vec + vec_read_in_block;
       bv.Close();
    }
    
    block_id = 0;
    vec_read_in_block = 0;
    queue->Complete();
    
    hts_close(out_vcf);
    bcf_hdr_destroy(out_hdr);
    out_open = false;
    std::cout << "BCF file with list of variant sites (" << arch_name + ".bcf" << ") created." << std::endl;
    
    
    if(bcf_index_build(out_vcf_name.c_str(), 14) == -1)
    {
        std::cout << "Error: index for BCF was not created." << std::endl;
        exit(1);
    }
     
    std::cout << "Index for BCF file with list of variant sites (" << arch_name + ".bcf.csi" << ") created." << std::endl;

    bcf_destroy1(rec);
    bcf_destroy1(new_rec);
    delete [] alt_desc;

    return true;
}


void VCFManager::addGTtoBitVector(bcf_hdr_t * hdr, bcf1_t * recc)
{
    int *gt_arr = NULL, ngt_arr = 0;
    int allele;
    bcf_get_genotypes(hdr, recc, &gt_arr, &ngt_arr);
    
    // Set vector with more significant bits of dibits
    for(int i = 0; i < ngt_arr; i++)
    {
        allele = bcf_gt_allele(gt_arr[i]);
        if(allele == 0 || allele == 1)
        {
            bv.PutBit(0);
        }
        else //if(bcf_gt_is_missing(gt_arr[i]) || bcf_gt_allele(gt_arr[i]) == 2)
        {
            bv.PutBit(1);
        }
    }
    bv.FlushPartialWordBuffer();
    
    // Set vector with less significant bits of dibits
    for(int i = 0; i < ngt_arr; i++)
    {
        allele = bcf_gt_allele(gt_arr[i]);
        if(allele == 1 || allele == 2)
        {
            bv.PutBit(1);
        }
        else //0
        {
            bv.PutBit(0);
        }
        
    }
    bv.FlushPartialWordBuffer();
    free(gt_arr);
    
    vec_read_in_block += 2; // Two vectors added
    if(vec_read_in_block == no_vec_in_block) // Insert complete block into queue of blocks
    {
        bv.TakeOwnership();
        
        queue->Push(block_id, bv.mem_buffer, vec_read_in_block);
        no_vec = no_vec + vec_read_in_block;
        block_id++;
        
        bv.Close();
        bv.Create(block_max_size);
        vec_read_in_block = 0;
    }
}
