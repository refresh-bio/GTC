/*
 This file is a part of GTC software distributed under GNU GPL 2 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */

#include "compression_settings.h"

char CompSettings::bits_used(unsigned int n) {
    char bits = 0;

    while(n)
    {
        n = n >> 1; bits++;
    }

    return bits;
}
