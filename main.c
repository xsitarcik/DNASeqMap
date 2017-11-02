#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwt.h"
#include "FMindex.h"
//#include "compression.h"

unsigned int genome_length;
unsigned int block_size = 128;
unsigned int sample_OCC_size = 128;
unsigned int sample_SA_size = 36;
unsigned char flag_compress = 1;
unsigned char flag_runs = 0;
unsigned char flag_mtf = 1;
unsigned char flag_huffman = 1;

int main ( int argc, char *argv[] )
{
 int i,j;
 unsigned int *suffix_array = NULL;
 unsigned int *sample_SA = NULL;
 unsigned char *s = "GGATCTGCGTCACAGCGAGCATAGCGAGAGCGGAGTTGCCGACGGCGGAGGCGACGCAGGGATCCGTCGGCCGCCATCCGCGGAAAGCATCCGCCCCCGAGGGGGACAGTCACTGACGCGGTCTTGCAGAGGCCTAGGGGGGCAGGTCAGTTTGAGTGGCTTGAGAACGCGAACTCTGGGATTACAGTGCAGTAATCTCCGGTCAACGGTGACGGCTTTAAGACAGGTCTTCGCAAAACCAGGCGGGGTGGCCTCGACGGGTTTTGCTGGTGGTTCAGGCGTACAATGCCCTGAAGAATAATTGAGAAAGTAGCACCCCTCGCCGCCTAGGATTACCTACCGGCGTCCGCCGCACCTTCGATTGTCGCGCCCACCCTCCCATTAGCCGGCACAGGTGGATGTGTCGCGACAGCCCGCCAAGATATCCTGAGGCGCAACGCGGACGGATGTCCCACGGAGTTGCCACAGGCGCCGAGCGCTTCACGGGCGACAGGAACTTG";
 unsigned char *bwt = NULL;
 unsigned int *bitvector_length;
 struct FMIndex *FM_index = NULL;
 struct compressedFMIndex *compressed_FM_index = NULL;
 unsigned char *alphabet = "ACGT";
 //load genome
 //abeceda
 //TODO

 genome_length = strlen(s);

 for (i=0;i<argc;i++)
  printf("%s\n",argv[i]);
 //count if memory is sufficient
 // for construction of SA we need cca 5n bytes 
 // then for bwt we need another n bytes
 // so together 6n bytes at least, 
 // if its not sufficient, break origin string TO DO

 //create suffix array and bwt of genome
 suffix_array = init_suffix_array(suffix_array,s, genome_length);
 bwt = create_bwt(suffix_array,s,genome_length);
 
 //free(s);
 
 //count if sample sizes are high enough so occ table can fit in memory TO DO
 //build FMIndex and free suffix_array AND GENOME (to do)
 if (flag_compress)
 {
  compressed_FM_index = build_compressed_FM_index(suffix_array,sample_SA_size, sample_OCC_size, genome_length,bwt,alphabet,flag_mtf, flag_runs, flag_huffman, block_size);

 }
 else
 {
  FM_index = build_FM_index(suffix_array,sample_SA_size, sample_OCC_size, genome_length,bwt,alphabet);
  
  printf("Org string is: %s\n",s);
  print_info_fm_index(FM_index);
///*
  int*result = search_pattern(FM_index,"AAAA");
  printf("results: %d %d\n",result[0],result[1]);
  while (result[0]<result[1])
  {
   printf("sa value %d je %d\n",result[0],get_SA_value(result[0],FM_index->bwt[result[0]],FM_index));
   result[0]++;
  }//*/
 }
 


 //app matching:
 //int*result = approximate_search(2,FM_index,"ACGAAACGATT");
 
 //align("TGTTAC","GGTTGAC", 2);


return 0;
}
