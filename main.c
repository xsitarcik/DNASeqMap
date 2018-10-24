#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bwt.h"
#include "FMindex.h"
#include <stdint.h>

//MAIN PARAMETERS:
//for constructing auxiliary tables of FM Index
unsigned int sample_OCC_size = 2; //in reality it's *64, MUST SET TO 2
unsigned int sample_SA_size = 32;
unsigned char max_error = 1;
unsigned int THRESHOLD = 500;
unsigned int k_mers_permutation = 9;
unsigned int total_kmers; //2^20 = 4^10
unsigned int *kmers_hash;
//>SRR493095.1 M00282:31:000000000-A0FFK:1:1:13945:1807 length=150



unsigned char load = 0;

//program parameters
unsigned char save = 0;
char *save_name;

char*filename_text;// = "alt_Celera_chr15.fa";
char *load_name;// = "alt_Celera_chr15_bwt_withoutN.txt";
char*filename_patterns;// = "SRR493095final.fasta";
//unsigned char*filename_patterns = "data/kratke_ready.fasta";

unsigned int MAX_READ_LENGTH = 200;
unsigned int READS_CHUNK = 70;

char *alphabet = "ACGT";
unsigned char alphabet_size = 3; //indexing from 0 so -1

unsigned char file_with_chunks = 1;
unsigned int CHUNK_SIZE = 70;


unsigned char flag_wavelet_tree = 0;
unsigned char flag_entries = 0;
unsigned char flag_use_gpu = 0;

//global program parametrrs
unsigned int genome_length;
 struct FMIndex_WT *FM_index_WT;
 struct wavelet_tree *WT_root;
unsigned char max_bits = sizeof(unsigned long long int)*8;
unsigned char*genome;
unsigned int*count_table;
unsigned int*entries;

int convert_string_to_int(char *string)
{
  int i,ret = 0;
  for(i=0; i<strlen(string); i++){
        ret = ret * 10 + ( string[i] - '0' );
      }
  return ret;
}

void print_help()
{ 
  printf("***Sequence Mapper - Master's thesis - by Jozef Sitarcik***\n");
  printf("Usage:\n");
  printf("SeqMap [Required input] [Required input parameter] [Optional input]*\n");
  printf(" Required input:\n");
  printf("  -r <string>   name of input file for Ref seqeunce\n");
  printf("  -p <string>   name of input file for Patterns\n");
  printf(" Required option - choose only one\n");
  printf("  -w            use FM Index based on Huffman-shaped Wavelet tree\n");
  printf("  -n            use FM Index aggregated in eNtries\n");
  printf("  -g            use GPU\n");
  printf(" Optional input:\n");
  printf("  -t <int>      max Threshold of values for aligning (default %d)\n",THRESHOLD);
  printf("  -e <int>      max allowed Error (default %d)\n", max_error);
  printf("  -i <string>   filename of prebuilt Index to load (defaultly not used)\n");
  printf("  -b <string>   filename to save Built Bwt (defaultly not used)\n");
  printf("  -c <int>      value for sampling Counters (default %d)\n",sample_OCC_size*64);
  printf("  -s <int>      value for sampling Suff arr values (default %d)\n",sample_SA_size);

}
int main ( int argc, char *argv[] )
{
 int i,j,k,count = 0;
 int ret;
 unsigned int *suffix_array = NULL;
 unsigned int *sample_SA = NULL;
 unsigned char *bwt = NULL;
 unsigned char c;
 unsigned int *bitvector_length;
 struct FMIndex *FM_index = NULL;
 unsigned int result_map;


 /*naive_compress("AAGCT");*/

//handle input options and parameters
for (i = 1; i<argc; i++)
{
  if (strcmp(argv[i],"-h") == 0)
  {
    print_help();
    return 0;
  }
  else if (strcmp(argv[i],"-t") == 0)
    THRESHOLD = convert_string_to_int(argv[++i]);
  else if (strcmp(argv[i],"-e") == 0)
    max_error = convert_string_to_int(argv[++i]);
  else if (strcmp(argv[i],"-s") == 0)
    sample_SA_size = convert_string_to_int(argv[++i]);
  else if (strcmp(argv[i],"-c") == 0){
    sample_OCC_size = convert_string_to_int(argv[++i])/64;
    if (sample_OCC_size<1){
      printf("Error - sampling size for counters cant be lower than 64 because of usage 64bitwords\n");
      return(-1);
    }
  }
  else if (strcmp(argv[i],"-w") == 0)
    flag_wavelet_tree = 1;
  else if (strcmp(argv[i],"-n") == 0)
    flag_entries = 1;
  else if (strcmp(argv[i],"-g") == 0){
    flag_use_gpu = 1;
    flag_entries = 1;
  }

  else if (strcmp(argv[i],"-r") == 0)
  {
    filename_text = realloc (filename_text, sizeof(char) * strlen(argv[++i])+1);
    strcpy(filename_text, argv[i]);
  }
  else if (strcmp(argv[i],"-p") == 0)
  {
    filename_patterns = realloc (filename_patterns,sizeof(char)*strlen(argv[++i])+1);
    strcpy(filename_patterns, argv[i]);
  }
  else if (strcmp(argv[i],"-i") == 0)
  {
    load = 1;
    load_name = realloc (load_name, sizeof(char)*strlen(argv[++i])+1);
    strcpy(load_name, argv[i]);
  }
  else if (strcmp(argv[i],"-b") == 0)
  {
    save = 1;
    save_name = realloc (save_name, sizeof(char)*strlen(argv[++i])+1);
    strcpy(save_name, argv[i]);
  }
}

 FILE *fp;
 FILE *fh_patterns;
 char *pattern = (char*)malloc(sizeof(char)*MAX_READ_LENGTH);
 

 if (load)
 {
  fp = fopen(load_name,"r");
  if (fp) {
   genome = load_genome_from_file_by_chunks(CHUNK_SIZE,filename_text,&genome_length);
   //genome = load_genome_from_file(filename_text,&genome_length);
   printf("Successfully loaded genome with length %u\n",genome_length);
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
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("construction of BWT took: %lf seconds\n", time_spent);
  //free(genome); 
  fflush(stdout);
 }
 
 fh_patterns = fopen(filename_patterns,"r");
 if (!(fh_patterns))
 {
  printf("Error when loading patterns from file %s\n",filename_patterns);
  exit(-1);
 }
 
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

if (flag_wavelet_tree){

  //build FM Index with WT for backward search
  FM_index_WT = build_FM_index_WT(suffix_array,bwt);
  
  unsigned int *result = (unsigned int*) malloc (sizeof(unsigned int)*2);
  if (result == NULL)
  {
   printf("error when allocating memory for result used when searching patterns\n");
   exit(-1);
  }

  printf("...approximate searching...\n");
  
  k = 0;
  count = 0;
  char o = 0;
  char line[256];
  clock_t begin1 = clock();
  fgets(line, 256, fh_patterns);
  //get name of the read
  while (1) 
  {
   pattern[0] = '\0';
   i = 0;
   while (1)
    {
     if (fgets(line, 256, fh_patterns)!= NULL);
      else{
        o = 1;
        break;
      }
     if (line[0] == '>')
      break;
     else {
      strncpy(&pattern[i],line, strlen(line)-1);
      i = i + strlen(line) - 1;
    }

    }
   if (o)
    break;

   /*printf("nacital: -%s-, %d\n",pattern,i);
   fflush(stdout);*/
   pattern[i-1]='\0';
   k++;
   if (strchr(pattern,'N')==NULL)
   {
    result_map = approximate_search_in_FM_index_WT(pattern,result);
    if (result_map)
    {
     ++count;
     printf("-%s-\n",pattern);
     fflush(stdout);
     printf("approximate position is %d\n",result_map);
    }
   }
  }
  clock_t end1 = clock();
  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("It took %lf seconds\n", time_spent1);
  printf("Total aligned reads: %d\n",count);
  printf("Total reads: %d\n",k);

  free(FM_index_WT->sampleSA);
  free(count_table);
  free(WT_root);
  free(FM_index_WT);
 }
 else if (flag_entries)
 {

  printf("...rebuilding FM Index into entries...\n");
  fflush(stdout);

  rebuild_FM_index_into_entries(suffix_array,bwt);

  unsigned int *result = (unsigned int*) malloc (sizeof(unsigned int)*2);
  if (result == NULL)
  {
   printf("error when allocating result array for searching patterns\n");
   exit(-1);
  }

  printf("...approximate searching...\n");
 k = 0;
  count = 0;
  char o = 0;
  char line[256];
  clock_t begin1 = clock();
  fgets(line, 256, fh_patterns);
  //get name of the read
  while (1) 
  {
   pattern[0] = '\0';
   i = 0;
   while (1)
    {
     if (fgets(line, 256, fh_patterns)!= NULL);
      else{
        o = 1;
        break;
      }
     if (line[0] == '>')
      break;
     else {
      strncpy(&pattern[i],line, strlen(line)-1);
      i = i + strlen(line) - 1;
    }

    }
   if (o)
   //if (k > 2500000)
    break;
   
   pattern[strlen(pattern)-1]='\0';
   k++;
   if (strchr(pattern,'N')==NULL)
   {
    result_map = approximate_search_in_FM_index_entry(pattern,result);
    if (result_map)
    {
     ++count; 
     printf("-%s-\n",pattern);
     fflush(stdout);
     printf("approximate position is %d\n",result_map);
    }
   }
  }
  //result_map = approximate_search_in_FM_index_entry("TCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCA",result);

  clock_t end1 = clock();
  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("It took %lf seconds\n", time_spent1);
  printf("Total aligned reads: %d\n",count);
  printf("Total reads: %d\n",k);
  
 }
 else 
 {
  printf("You need to specify -g -w or -n\n");
  exit(-1);  
  }
return 0;
}
