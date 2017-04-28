#!/bin/bash
# This script removes all dependencies of GTC

rm -fr sdsl-lite*
rm -fr htslib-1.3.2*
rm -fr share
rm -fr lib
rm include/divsufsort*
rm -rf  include/htslib
rm -rf  include/sdsl

rm bin/bgzip
rm bin/htsfile
rm bin/tabix

rm -rf bin
