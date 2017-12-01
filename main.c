#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bwt.h"
#include "FMindex.h"
#include <stdint.h>

#include "compression.h"

unsigned int genome_length;
unsigned int block_size = 4000;
unsigned int sample_OCC_size = 64;
unsigned int sample_SA_size = 32;
unsigned char flag_compress = 1;
unsigned char flag_runs = 3;
unsigned char flag_mtf = 1;
unsigned char flag_huffman = 1;
unsigned char max_errors=0;


int main ( int argc, char *argv[] )
{
 int i,j;
 unsigned int *suffix_array = NULL;
 unsigned int *sample_SA = NULL;
 unsigned char *genome; 
 unsigned char *bwt = NULL;
 unsigned int *bitvector_length;
 struct FMIndex *FM_index = NULL;
 struct compressedFMIndex *compressed_FM_index = NULL;
 //unsigned char *alphabet = "ACGT";
 unsigned char *alphabet = "ACGT";
 unsigned char*filename = "testik.txt";
 //load genome
 //abeceda
 //TODO
 
 for (i=0;i<argc;i++)
 {
  printf("%s\n",argv[i]);
 }

 genome = load_genome_from_file(filename,&genome_length);
 //reverse_string(genome);
 //count if memory is sufficient
 // for construction of SA we need cca 5n bytes 
 // then for bwt we need another n bytes
 // so together 6n bytes at least, 
 // if its not sufficient, break origin string TO DO

 //create suffix array and bwt of genome
 clock_t begin = clock();
 suffix_array = init_suffix_array(suffix_array,genome, genome_length);
 bwt = create_bwt(suffix_array,genome,genome_length);
  clock_t end = clock();
 double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
 printf("trvanie kompresie: %lf", time_spent);
 //free(s);

 //count if sample sizes are high enough so occ table can fit in memory TO DO
 //build FMIndex and free suffix_array AND GENOME (to do)
 if (flag_compress)
 {
  compressed_FM_index = build_compressed_FM_index(suffix_array,sample_SA_size, sample_OCC_size, genome_length,bwt,alphabet,flag_mtf, flag_runs, flag_huffman, block_size);
  
  /*int*result = approximate_search_in_compressed_FM_index(max_errors,compressed_FM_index,"TATTATATTAATTATCATCCTAACTGAGG",flag_mtf,flag_runs,flag_huffman);
  printf("results: %d %d\n",result[0],result[1]);
  printf("ide sa na vypocet\n");
  fflush(stdout);
  printf("wata %d",__builtin_popcount(500));*/

  for (i=0;i<compressed_FM_index->length;i++)
  {
   free(compressed_FM_index->array_of_blocks[i].occurences);
   free(compressed_FM_index->array_of_blocks[i].bitvector);
   free_huffman_tree(compressed_FM_index->array_of_blocks[i].huffman_tree);
  }
  free(compressed_FM_index->array_of_blocks);
  free(compressed_FM_index->count_table);
  free(compressed_FM_index->sampleSA);
  free(compressed_FM_index);
 }
 else
 {
  FM_index = build_FM_index(suffix_array,sample_SA_size, sample_OCC_size, genome_length,bwt,alphabet);
  
  printf("Org string is: %s\n",genome);
  print_info_fm_index(FM_index);
///*

  int*result = approximate_search(max_errors,FM_index,"TYYASI");
  printf("results: %d %d\n",result[0],result[1]);
  while (result[0]<result[1])
  {
   printf("sa value %d je %d\n",result[0],get_SA_value(result[0],FM_index->bwt[result[0]],FM_index));
   result[0]++;
  }//*/


   free(FM_index->sampleSA);
   free(FM_index->count_table);
   for (i=0;i<strlen(alphabet);i++)
    free(FM_index->occurence_table[i]);
   free(FM_index->occurence_table);
   free(FM_index->bwt);
   free(FM_index);
 }
 


 //app matching:
 //int*result = approximate_search(2,FM_index,"ACGAAACGATT");
 
 //align("TGTTAC","GGTTGAC", 2);



return 0;
}
