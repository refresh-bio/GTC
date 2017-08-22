# GTC - GenoTypes Compressor

Requirements
--------------

GTC requires:

* A modern, C++11 ready compiler such as `g++` version 4.9 or higher or `clang` version 3.2 or higher.
* The CMake build system.
* A 64-bit operating system. Either Mac OS X or Linux are currently supported.
* For best performance the processor of the system should support fast bit operations available in `SSE4.2`


Installation
--------------

To download, build and install GTC use the following commands.
```sh
git clone https://github.com/refresh-bio/GTC.git
cd GTC
./install.sh 
make
make install
```
The `install.sh` script downloads and installs the SDSL and HTSlib libraries into the `include` and `lib` directories in the `GTC` directory. 

By default GTC is installed in the `bin` directory of the `/usr/local/` directory. A different location prefix can be specified with `prefix` parameter:
```sh
make prefix=/usr/local install
```
---
To uninstall GTC:
```sh
make uninstall
```
This uninstalls GTC from the `/usr/local` directory. To uninstall from different location use the `prefix` parameter:
```sh
make prefix=/usr/local uninstall
```
To uninstall the SDSL and HTSlib libraries use the provided uninstall script:
```sh
./uninstall.sh 
```
---
To clean the GTC build use:
```sh
make clean
```
Usage
--------------
* Compress the input VCF/BCF file
```
Input: [file_name] VCF/BCF file. 
Output: [archive_name].ind file with samples names, [archive_name].bcf file with variant sites description, [archive_name].gtc file with the archive. By default [archive_name] is set to "archive".

Usage: gtc compress <options> [file_name] 
[file_name]		- input file (a VCF or VCF.GZ file by default)

Available options (optional): 
Input: 
	-b    	- input is a BCF file (input is a VCF or VCF.GZ file by default)	
	-p [x]	- set ploidy of samples in input VCF to [x] (number >= 1; 2 by default)
Output: 
	-o [name]	- set archive name to [name]("archive" by default)	
Parameters: 
	-t [x]	- set number of threads to [x] (number >= 1; 8 by default)
	-d [x]	- set maximum depth to [x] (number >= 0; 0 means no matches; 100 by default)
	-g [x]   	- [DEV] set number of vector groups [percentage of 1s] to [x] (max: 32; 32 by default)	
	-hm [x]   	- [DEV] set n_vec_history for matches to pow(2, [x]) (9 by default, min: 8)	
	-hc [x]   	- [DEV] set n_vec_history for copies to pow(2, [x]) (17 by default, min: 8)	
  ```
  
 * Decompress / Query the archive.
 ```
Input [archive_name] archive ([archive_name].ind, [archive_name].bcf and [archive_name].gtc). 
Output: a VCF/BCF file.

Usage: gtc view <options> [archive_name]
Available options: 
Output: 
	-o [name]	- output to a file and set output name to [name] (stdout by default)	
	-b	- output a BCF file (output is a VCF file by default)	
	-C 	- write AC/AN to the INFO field (always set when using -minAC, -maxAC, -minAF or -maxAF)
	-G 	- don't output sample genotypes (only #CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO columns)
	-c [0-9]   set level of compression of the output bcf (number from 0 to 9; 1 by default; 0 means no compression)	
Query: 
	-r	- range in format [chr]:[start]-[end] (for example: -r 14:19000000-19500000). By default all variants are decompressed.
	-s	- sample name(s), separated by comms (for example: -s HG00096,HG00097) OR '@' sign followed by the name of a file with sample name(s) separated by whitespaces (for exaple: -s @file_with_IDs.txt). By default all samples/individuals are decompressed
	-n X 	- process at most X records (by default: all from the selected range)
Settings: 
	-minAC X 	- report only sites with count of alternate alleles among selected samples smaller than or equal to X (default: no limit)
	-maxAC X 	- report only sites with count of alternate alleles among selected samples greater than or equal to X
	-minAF X 	- report only sites with allele frequency among selected samples greather than or equal to X (X - number between 0 and 1; default: 0)
	-maxAF X 	- report only sites with allele frequency among selected samples smaller than or equal to X (X - number between 0 and 1; default: 1)
	-m X	- limit maximum memory usage to remember previous vectors to X MB (no limit by default)	
 ```

Toy example
--------------

There is an example VCF file, `toy.vcf`, in the `toy_ex` folder, wchich can be used to test GTC

To compress the example VCF file and store the archive called `toy_arch` in the `toy_ex` folder:
```sh
./gtc compress -o toy_ex/toy_arch toy_ex/toy.vcf
```
This will create an archive consisting of four files:
* `toy_arch.bcf` - BCF file with all variant sites description,
* `toy_arch.bcf.csi` - BCF index file,
* `toy_arch.gtc` - main archive of all genotypes,
* `toy_arch.ind` - list of all individuals.

To view the compressed archive (to decompress it) in VCF format:
```sh
./gtc view toy_ex/toy_arch
```

For more options see Usage section.


Developers
--------------
The GTC algorithm was invented by [Agnieszka Danek](https://github.com/agnieszkadanek) and Sebastian Deorowicz.
The implementation is by Agnieszka Danek (mainly) and Sebastian Deorowicz.

Publication
--------------
Danek, A., Deorowicz, S., GTC: an attempt to maintenance of huge genome collections compressed, [bioRxiv](http://biorxiv.org/content/early/2017/04/28/131649), 2017;

