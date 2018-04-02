#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bwt.h"
#include "FMindex.h"
#include <stdint.h>
#include "compression.h"

unsigned int genome_length;
unsigned char max_error = 1;
unsigned char max_bits = sizeof(unsigned long long int)*8;
//MAIN PARAMETERS:
//for constructing auxiliary tables of FM Index
unsigned int sample_OCC_size = 80; //in reality it's *64
unsigned int sample_SA_size = 32;

//program parameters
unsigned char save = 0;
unsigned char *save_name = "alt_Celera_chr15_bwt.txt";
unsigned char load = 1;
unsigned char *load_name = "alt_Celera_chr15_bwt.txt";

unsigned char*filename_text = "alt_Celera_chr15.fa";
unsigned char*filename_patterns = "e_coli_10000snp.fa";
unsigned int MAX_READ_LENGTH = 60;
unsigned char *alphabet = "ACGNT";
unsigned char alphabet_size = 4; //indexing from 0

unsigned char file_with_chunks = 1;
unsigned int CHUNK_SIZE = 70;

//for compression
unsigned int block_size = 15000;
unsigned char flag_compress = 0;
unsigned char flag_runs = 7;
unsigned char flag_mtf = 1;
unsigned char flag_huffman = 1;
unsigned char flag_wavelet_tree = 1;

//for string matching


int main ( int argc, char *argv[] )
{
 int i,j,k,count = 0;
 unsigned int *suffix_array = NULL;
 unsigned int *sample_SA = NULL;
 unsigned char *genome; 
 unsigned char *bwt = NULL;
 unsigned char c;
 unsigned int *bitvector_length;
 struct FMIndex *FM_index = NULL;
 struct compressedFMIndex *compressed_FM_index = NULL;
 struct FMIndex_WT *FM_index_WT = NULL;

 FILE *fp;
 FILE *fh_patterns;
 char *pattern = (char*)malloc(sizeof(char)*MAX_READ_LENGTH);
 /*
 for (i=0;i<argc;i++)
 {
  printf("%s\n",argv[i]);
 }*/

 if (load)
 {
  fp = fopen(load_name,"r");
  if (fp) {
   printf("...loading BWT from file %s... \n",load_name);
   genome_length = 0;
   fscanf (fp, "%u\n", &genome_length); 
   printf("genome length : %d\n",genome_length);
   bwt = (unsigned char*)malloc(sizeof(unsigned char)*genome_length+1);
   fread (bwt, 1, genome_length, fp);
   bwt[genome_length]='\0';
    
   getc(fp);//read newline
   suffix_array = (unsigned int*)malloc(sizeof(unsigned int)*genome_length);
   for (i=0;i<genome_length;i++)
    fscanf (fp, "%u,", &suffix_array[i]);
   fclose(fp);
  }
 }
 else 
 {
  //load main string from file
  if (file_with_chunks){
   genome = load_genome_from_file_by_chunks(CHUNK_SIZE,filename_text,&genome_length);
   printf("loading by chunks\n");
  }
  else{
   genome = load_genome_from_file(filename_text,&genome_length);
   printf("normal loading\n");
  }
  if (genome_length<=1)
  {
   printf("Error when reading file: %s\n",filename_text);
   exit(-1);
  }
  else
  printf("Size of read genome is %d characters\n",genome_length);

  //printf("read%s\n",genome);
  //create suffix array and bwt of main input string
  clock_t begin = clock();
  printf("...constructing BWT...\n");
  suffix_array = init_suffix_array(suffix_array,genome,genome_length);
  
  bwt = create_bwt(suffix_array,genome,genome_length);
  //printf("%s\n",bwt);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("construction of BWT took: %lf seconds\n", time_spent);
  free(genome); 
  fflush(stdout);
 }
 
 fh_patterns = fopen(filename_patterns,"r");
 if (!(fh_patterns))
 {
  printf("Error when loading patterns from file %s\n",filename_patterns);
  exit(-1);
 }

 //reverse_string(genome);

 //should handle cases like:
 //count if memory is sufficient
 // for construction of SA we need cca 5n bytes 
 // then for bwt we need another n bytes
 // so together 6n bytes at least, 
 // if its not sufficient, break origin string

 
 
//save created BWT in file
if (save)
{
 fp = fopen(save_name,"w");
 if(fp == NULL)
 {
  printf("\nProblem when creating file %s for saving index!",save_name);
  exit(-1);
 }
 fprintf(fp,"%u\n",genome_length);
 fprintf(fp,"%s\n",bwt);
 //fprintf(fp,"");
 for (i=0;i<genome_length;i++)
  fprintf(fp,"%u,",suffix_array[i]);
    
 fclose(fp);
 }


 //build FMIndex and free suffix_array AND GENOME (to do)
 if (flag_compress)
 {
  compressed_FM_index = build_compressed_FM_index(suffix_array,bwt,flag_mtf, flag_runs, flag_huffman, block_size);
   printf("ide sa na vypocet\n");
   fflush(stdout);
  int*result = approximate_search_in_compressed_FM_index(compressed_FM_index,"TAATCGGTGGGAGTATTCAACGTGATGAAGAC",flag_mtf,flag_runs,flag_huffman);
  printf("results: %d %d\n",result[0],result[1]);
  printf("ide sa na vypocet\n");
  fflush(stdout);
  

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

  //build FM Index with WT for backward search
  FM_index_WT = build_FM_index_WT(suffix_array,bwt);
 
  printf("--------------------------------------\n");
  
  unsigned int*result = (unsigned int*)malloc(2*sizeof(unsigned int));
  printf("...approximate searching...\n");
  k = 0;
  clock_t begin1 = clock();
  while (fgets(pattern, MAX_READ_LENGTH, fh_patterns) != NULL) {
    fgets(pattern, MAX_READ_LENGTH, fh_patterns); //in case 
    pattern[strlen(pattern)-1]='\0';
    /*printf("--%s--\n", pattern);
    fflush(stdout);*/
    k++;
    i = 0; j = 0;
    while (i!=strlen(pattern))
      if (pattern[i++]=='N')
        j++;
    if (j<=max_error)
     count += approximate_search_in_FM_index_WT(FM_index_WT,pattern,result);
    //else
      //printf("too much Ns\n");
  }
  clock_t end1 = clock();
  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("It took %lf seconds\n", time_spent1);
  printf("Total aligned reads: %d\n",count);
  printf("Total reads: %d\n",k);
  //long long int result = approximate_search_in_FM_index_WT(max_errors,FM_index_WT,"TCGATTATATCACTTAATGACTTTTGGGTCAGGGTGTGTTACCTTACAGGAATTGAGACCGTCCATTAATTTCTCTTGCATTTAT");
  //printf("result je %lld\n",result);
  //build FM Index with WT for forward search
   //..
//CGTAACTCG-TT GCTCAGTCGCTTT TCGTACTGCGCG -2errors
//CGTATCTCGATT GCTCAGTCGCTTT TCGTACTGCGCG withour error
  /*
  
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
  FM_index = build_FM_index(suffix_array,bwt);
  
  print_info_fm_index(FM_index);
  int*result = approximate_search(FM_index,"TYYASI");
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
