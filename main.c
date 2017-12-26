#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bwt.h"
#include "FMindex.h"
#include <stdint.h>
#include "compression.h"

unsigned int genome_length;

//MAIN PARAMETERS:
//for constructing auxiliary tables of FM Index
unsigned int sample_OCC_size = 256;
unsigned int sample_SA_size = 32;

//for compression
unsigned int block_size = 1000;
unsigned char flag_compress = 0;
unsigned char flag_runs = 7;
unsigned char flag_mtf = 1;
unsigned char flag_huffman = 1;
unsigned char flag_wavelet_tree = 1;

//for string matching
unsigned char max_errors=0;

int main ( int argc, char *argv[] )
{
 int i,j,k;
 unsigned int *suffix_array = NULL;
 unsigned int *sample_SA = NULL;
 unsigned char *genome; 
 unsigned char *bwt = NULL;
 unsigned int *bitvector_length;
 struct FMIndex *FM_index = NULL;
 struct compressedFMIndex *compressed_FM_index = NULL;
 struct FMIndex_WT *FM_index_WT = NULL;
 unsigned char *alphabet = "ACGT";
 unsigned char*filename = "test.txt";
 
 /*
 for (i=0;i<argc;i++)
 {
  printf("%s\n",argv[i]);
 }*/
 
 //load main string from file
 genome = load_genome_from_file(filename,&genome_length);
 
 //reverse_string(genome);

 //should handle cases like:
 //count if memory is sufficient
 // for construction of SA we need cca 5n bytes 
 // then for bwt we need another n bytes
 // so together 6n bytes at least, 
 // if its not sufficient, break origin string

 //create suffix array and bwt of main input string
clock_t begin = clock();
printf("...constructing BWT...");
suffix_array = init_suffix_array(suffix_array,genome, genome_length);
bwt = create_bwt(suffix_array,genome,genome_length);
clock_t end = clock();
double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
printf("construction of BWT took: %lf seconds", time_spent);
free(genome); 
fflush(stdout);
 /*
   //save created BWT in file
    FILE *fp;
    fp = fopen("bwt.txt","w");
    if(fp == NULL)
    {
        printf("\nSorry! Couldn't Create File");
        exit(-1);
    }
    fprintf(fp,"%s",bwt);
    fclose(fp);
    */


 //build FMIndex and free suffix_array AND GENOME (to do)
 if (flag_compress)
 {
  compressed_FM_index = build_compressed_FM_index(suffix_array,sample_SA_size, sample_OCC_size, genome_length,bwt,alphabet,flag_mtf, flag_runs, flag_huffman, block_size);
  
  /*int*result = approximate_search_in_compressed_FM_index(max_errors,compressed_FM_index,"TATTATATTAATTATCATCCTAACTGAGG",flag_mtf,flag_runs,flag_huffman);
  printf("results: %d %d\n",result[0],result[1]);
  printf("ide sa na vypocet\n");
  fflush(stdout);
  printf("wata %d",__builtin_popcount(500));*/
  

  /*clock_t begin2 = clock();
  for (i=0;i<10000;i=i+10){
    j = i/block_size;
    k = i - j*block_size;
    unsigned char *ret = decompress_block(compressed_FM_index->array_of_blocks[j].bitvector_length,compressed_FM_index->array_of_blocks[j].bitvector,flag_mtf,flag_runs,flag_huffman,k+1,alphabet,compressed_FM_index->array_of_blocks[j].huffman_tree);
    free(ret);
  }
   clock_t end2 = clock();
  double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
  printf("%d access operacie %lf", i/10, time_spent2);
  */

  //freeing FM Index
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
 else if (flag_wavelet_tree){
  FM_index_WT = build_FM_index_WT(suffix_array,sample_SA_size,sample_OCC_size,genome_length,bwt,alphabet);
  /*
  clock_t begin1 = clock();
  for (i=0;i<genome_length;i++){
   wt_access(i,FM_index_WT->WT_root,sample_OCC_size);
   //printf("na pozicii %d je znak %c\n",i,wt_access(i,FM_index_WT->WT_root,sample_OCC_size));
   //printf("pocet %c do pozicie %d je %d\n",'C',69,wt_rank('A',69+1,root,sample_OCC_size));
  }
  clock_t end1 = clock();
  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("%d access operacie %lf", i, time_spent1);
  */
 }
 else
 {
  FM_index = build_FM_index(suffix_array,sample_SA_size, sample_OCC_size, genome_length,bwt,alphabet);
  
  print_info_fm_index(FM_index);
  int*result = approximate_search(max_errors,FM_index,"TYYASI");
  printf("results: %d %d\n",result[0],result[1]);
  while (result[0]<result[1])
  {
   printf("sa value %d je %d\n",result[0],get_SA_value(result[0],FM_index->bwt[result[0]],FM_index));
   result[0]++;
  }
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
