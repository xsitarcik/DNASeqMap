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
unsigned int k_mers_permutation = 1;
unsigned int total_kmers; //2^20 = 4^10
unsigned int *kmers_hash;

//program parameters
unsigned char load = 0;
unsigned char save = 0;
char *save_name;
char*filename_text;
char *load_name;
char*filename_patterns;

unsigned int MAX_READ_LENGTH = 250;

char *alphabet = "ACGT";
unsigned char alphabet_size = 3; //indexing from 0 so -1

unsigned char flag_wavelet_tree = 0;
unsigned char flag_entries = 0;
unsigned char flag_rankbitvectors=0;
unsigned char flag_withoutSA = 0;
unsigned char flag_locate = 0;
unsigned char flag_count = 0;

//global program parametrrs
unsigned int genome_length;
 struct FMIndex_WT *FM_index_WT;
 struct wavelet_tree *WT_root;
unsigned char max_bits = sizeof(unsigned long long int)*8;
unsigned char*genome;
unsigned int*count_table;
unsigned int*entries;

unsigned int*root_rankvector;
unsigned int*left_rankvector;
unsigned int*right_rankvector;
unsigned int *sample_SA;

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
  printf("***FM Index improvements for DNA sequences - by Jozef Sitarcik***\n");
  printf("Usage:\n");
  printf("SeqMap [Required input] [Required input parameter] [Optional input]*\n");
  printf(" Required input:\n");
  printf("  -r <string>   name of input file for Ref seqeunce\n");
  printf("  -p <string>   name of input file for Patterns\n");
  printf(" Required option - choose only one\n");
  printf("  -w            use FM Index based on Huffman-shaped Wavelet tree\n");
  printf("  -v            use FM Index based on bitvector interleaved with ranks\n");
  printf("  -n            use FM Index aggregated in eNtries\n");
  printf("  -g            use FM Index aggregated in eNtries whitout suffix array\n");
  printf(" Optional input:\n");
  printf("  -t <int>      max Threshold of values for aligning (default %d)\n",THRESHOLD);
  printf("  -e <int>      max allowed Error (default %d)\n", max_error);
  printf("  -i <string>   filename of prebuilt Index to load (defaultly not used)\n");
  printf("  -b <string>   filename to save Built Bwt (defaultly not used)\n");
  printf("  -c <int>      value for sampling Counters (default %d)\n",sample_OCC_size*64);
  printf("  -s <int>      value for sampling Suff arr values (default %d)\n",sample_SA_size);
  printf("  -k <int>      value of k in kmer hash table (default %d)\n",k_mers_permutation);
  printf("  -x            run only operation locate for input(defaultly not used)\n");
  printf("  -y            run only operation count for input(defaultly not used)\n");

}
int main ( int argc, char *argv[] )
{
 int i,j,k,count = 0;
 int ret;
 unsigned int *suffix_array = NULL;
 unsigned char *bwt = NULL;
 unsigned char c;
 unsigned int *bitvector_length;
 struct FMIndex *FM_index = NULL;
 unsigned int result_map;

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
  else if (strcmp(argv[i],"-k") == 0)
    k_mers_permutation = convert_string_to_int(argv[++i]);
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
  else if (strcmp(argv[i],"-v") == 0)
    flag_rankbitvectors = 1;
  else if (strcmp(argv[i],"-g") == 0)
    flag_withoutSA = 1;
  else if (strcmp(argv[i],"-x") == 0)
    flag_locate = 1;
  else if (strcmp(argv[i],"-y") == 0)
    flag_count = 1;

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
   genome = load_ref_sequence(filename_text,&genome_length);
   printf("...loading BWT from file %s... \n",load_name);
   genome_length = 0;
   fscanf (fp, "%u\n", &genome_length); 
   printf("genome length : %d\n",genome_length);
   bwt = (unsigned char*)malloc(sizeof(unsigned char)*genome_length+1);
   fread (bwt, 1, genome_length, fp);
   bwt[genome_length]='\0';
    
   getc(fp);//read newline
   //suffix_array = (unsigned int*)malloc(sizeof(unsigned int)*genome_length);
   unsigned int samples_count = (genome_length+sample_SA_size)/sample_SA_size;
   sample_SA = (unsigned int *) malloc((samples_count)*sizeof(unsigned int));
   for (i=0;i<samples_count;i++)
    fscanf (fp, "%u,", &sample_SA[i]);
   fclose(fp);
   printf("nacitanych bolo %u\n",samples_count);
  }
 }
 else 
 {
  //load main string from file
  genome = load_ref_sequence(filename_text,&genome_length);
  if (genome_length<=1)
  {
   printf("Error when reading file: %s\n",filename_text);
   exit(-1);
  }
  else
  printf("Reference sequence contains %d characters\n",genome_length);

  //create suffix array and bwt of main input string
  clock_t begin = clock();
  printf("...constructing BWT...\n");
  suffix_array = init_suffix_array(suffix_array,genome,genome_length);
  bwt = create_bwt(suffix_array,genome,genome_length);
  sample_SA = create_sample_SA(suffix_array);
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("construction of BWT took: %lf seconds\n", time_spent);
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
 for (i=0;i<genome_length;i=i+sample_SA_size)
  fprintf(fp,"%u,",sample_SA[i/64]);
 printf("ulozenych bolo %u\n",i);   
 fclose(fp);
 }

if (flag_wavelet_tree){

  //build FM Index with WT for backward search
  FM_index_WT = build_FM_index_WT(suffix_array,bwt);
  unsigned int *result = (unsigned int*) malloc (sizeof(unsigned int)*2);
  if (result == NULL){
   printf("error when allocating memory for result used when searching patterns\n");
   exit(-1);
  }

  printf("...approximate searching...\n");
  k = 0; count = 0;
  char o = 0; char line[256];
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

   pattern[i-1]='\0';
   k++;
   if (strchr(pattern,'N')==NULL)
   {
    result_map = approximate_search_in_FM_index_WT(pattern,result);
    if (result_map)
    {
     ++count;
     printf("-%s-\n",pattern);
     printf("approximate position is %d\n",result_map);
    }
   }
  }
  clock_t end1 = clock();
  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("It took %lf seconds\n", time_spent1);
  printf("Total aligned reads: %d\n",count);
  printf("Total reads: %d\n",k+1);

  free(FM_index_WT->sampleSA);
  free(count_table);
  free(WT_root);
  free(FM_index_WT);
 }
 else if (flag_entries)
 {

  printf("...rebuilding FM Index into entries...\n");
  fflush(stdout);

  rebuild_FM_index_into_entries(suffix_array,bwt,flag_withoutSA);

  unsigned int *result = (unsigned int*) malloc (sizeof(unsigned int)*2);
  if (result == NULL)
  {
   printf("error when allocating result array for searching patterns\n");
   exit(-1);
  }

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
   
   pattern[strlen(pattern)]='\0';
   k++;

   if (flag_count)
    result_map = search_pattern_in_FM_index_entry(pattern,strlen(pattern),result);  
   else if(flag_locate){
    result_map = search_pattern_in_FM_index_entry(pattern,strlen(pattern),result);  
    for (i=result[0];i<result[1];i++){
      get_SA_value_entry(i);
      ++count;
      }
    }
   else{
    result_map = approximate_search_in_FM_index_entry(pattern,result);
    if (result_map)
    {
     printf("%s\n",pattern);
     printf("approximate position is %u\n",result_map);
     ++count; 
    }
   }
  }

  clock_t end1 = clock();
  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("It took %lf seconds\n", time_spent1);
  printf("Total aligned reads: %d\n",count);
  printf("Total reads: %d\n",k+1);
 }
 
 else if (flag_rankbitvectors){

  printf("...rankbitvector FM Index...\n");
  fflush(stdout);

  create_rank_bitvectors(bwt);
  unsigned int *result = (unsigned int*) malloc (sizeof(unsigned int)*2);

  k = 0;
  unsigned long counter2 = 0;
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
   
   pattern[strlen(pattern)]='\0';
   k++;

   if (flag_count)
    result_map = search_pattern_in_rankbitvector(pattern,strlen(pattern),result);  
   else if(flag_locate){
    result_map = search_pattern_in_rankbitvector(pattern,strlen(pattern),result);  
    for (i=result[0];i<result[1];i++){
      get_SA_value_rankbitvector(i);
      ++counter2;
      }
    }
   else{
    /*result_map = approximate_search_in_FM_index_entry(pattern,result);
    if (result_map)
    {
     printf("%s\n",pattern);
     printf("approximate position is %u\n",result_map);
     ++count; 
    }*/
    printf("error. Approximate search is not yet supported for rankbitvectors\n");
    exit(-1);
   }
  }
 
  clock_t end1 = clock();
  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("It took %lf seconds\n", time_spent1);
  printf("Total reads: %d\n",k+1);
  printf("Total occurences processed: %lu\n",counter2);
  printf("You need to specify -g -w or -n\n");
  exit(-1);  
  }
  else if (flag_withoutSA){

    printf("...rebuilding FM Index into entries without SA\n");
    rebuild_FM_index_into_entries(suffix_array,bwt,flag_withoutSA);

    unsigned int *result = (unsigned int*) malloc (sizeof(unsigned int)*2);
    if (result == NULL)
    {
      printf("error when allocating result array for searching patterns\n");
      exit(-1);
    }

    k = 0;  count = 0;
    char o = 0; char line[256];
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
   
   pattern[strlen(pattern)]='\0';
   k++;

   if (flag_count)
    result_map = search_pattern_in_FM_index_exSA(pattern,strlen(pattern),result);  
   else if(flag_locate){
    result_map = search_pattern_in_FM_index_exSA(pattern,strlen(pattern),result);  
    for (i=result[0];i<result[1];i++){
      get_SA_value_exSA(i);
      ++count;
      }
    }
   else{
    /*result_map = approximate_search_in_FM_index_entry(pattern,result);
    if (result_map)
    {
     printf("%s\n",pattern);
     printf("approximate position is %u\n",result_map);
     ++count; 
    }*/
    printf("error. Approximate search is not yet supported this option\n");
    exit(-1);
   }

  }

  clock_t end1 = clock();
  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("It took %lf seconds\n", time_spent1);
  printf("Total aligned reads: %d\n",count);
  printf("Total reads: %d\n",k+1);
 }


  else {
    printf("Error. Make sure that you selected mandatory option for a type of FM Index. Use -h to show additional help.");
    exit(-1);
  }

return 0;
}
