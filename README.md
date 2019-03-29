# DNASeqMap
This is a cache-efficient implementation of FM-index focused solely on DNA sequences. Program allows to approximately search input patterns in a reference sequence.

## Input parameters
`DNASeqMap [Required input] [Required input parameter] [Optional input]*`   
Required input:   
&nbsp;&nbsp;  `-r <string>   name of input file for Ref seqeunce`   
&nbsp;&nbsp;  `-p <string>   name of input file for Patterns`   
Required option - choose only one    
&nbsp;&nbsp;  `-v            use FM Index based on WT bitvectors interleaved with ranks`   
&nbsp;&nbsp;  `-n            use FM Index interleaving all the data in chunks`   
&nbsp;&nbsp;  `-g            use FM Index interleaved all the data except suffix array`   
Optional input:    
&nbsp;&nbsp;  `-t <int>      max Threshold of values for aligning (default 500)`   
&nbsp;&nbsp;  `-e <int>      max allowed Error (default 1)`   
&nbsp;&nbsp;  `-i <string>   filename of prebuilt Index to load (defaultly not used)`   
&nbsp;&nbsp;  `-b <string>   filename to save Built Bwt (defaultly not used)`      
&nbsp;&nbsp;  `-k <int>      value of k in kmer hash table (default 1)`   
&nbsp;&nbsp;  `-x            run only operation locate for input(defaultly not used)`   
&nbsp;&nbsp;  `-y            run only operation count for input(defaultly not used)`   

## Compilation
DNASeqMap file is already executable compiled file. Otherwise compilation is as follows (with use of [Clang compiler](https://clang.llvm.org/)):   
`clang main.c bwt.c FMindex.c wavelet_tree.c -o DNASeqMap -O3 -funroll-loops -fomit-frame-pointer -ffast-math`

## Examples
To map input sequences from file `patterns.fa` in reference sequence from file `refsequence.fa` and use FM-index interleaved without sampled suffix array:   
&nbsp;&nbsp;`./DNASeqMap -r refsequence.fa -p patterns.fa -g`

To run only operation **locate(P)** use `-x`, e.g.:  
&nbsp;&nbsp;`./DNASeqMap -r refsequence.fa -p patterns.fa -g -x`   

To run only operation **count(P)** use `-y`, e.g.:  
&nbsp;&nbsp;`./DNASeqMap -r refsequence.fa -p patterns.fa -g -y`  

To use k-mer table, it is needed to set k to some value, e.g. to set k=7:   
&nbsp;&nbsp;`./DNASeqMap -r refsequence.fa -p patterns.fa -g -k 7`   
