#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bwt.h"
#include "FMindex.h"
#include <stdint.h>
#include "compression.h"
#include <CL/cl.h>

#define DEVICE_GROUP_SIZE 256
#define NUM_OF_READS_PER_WARP 8
#define WARP_COUNT 3
#define MAX_SOURCE_SIZE 0x100000
//MAIN PARAMETERS:
//for constructing auxiliary tables of FM Index
unsigned int sample_OCC_size = 2; //in reality it's *64, MUST SET TO 2
unsigned int sample_SA_size = 32;
unsigned char max_error = 1;
int MAX_RESULTS = 200;
int MIN_RESULT_LENGTH = 40; //pattern length /(maxerrror+1) / 2;
unsigned int THRESHOLD = 150;

//>SRR493095.1 M00282:31:000000000-A0FFK:1:1:13945:1807 length=150

//program parameters
unsigned char save = 0;
unsigned char *save_name = "alt_Celera_chr15_bwt_withoutN.txt";

unsigned char load = 1;

unsigned char *load_name = "patdesiattisic_bwt2.txt";
unsigned char*filename_text = "patdesiattisic.txt";
unsigned char*filename_patterns = "meko.fa";

/*unsigned char*filename_text = "alt_Celera_chr15.fa";
unsigned char *load_name = "alt_Celera_chr15_bwt_withoutN.txt";
unsigned char*filename_patterns = "SRR493095final.fasta";
*/

unsigned int MAX_READ_LENGTH = 200;
unsigned int READS_CHUNK = 70;

char *alphabet = "ACGT";
unsigned char alphabet_size = 3; //indexing from 0 so -1

unsigned char file_with_chunks = 1;
unsigned int CHUNK_SIZE = 70;

//for compression
unsigned int block_size = 15000;
unsigned char flag_compress = 0;
unsigned char flag_runs = 7;
unsigned char flag_mtf = 1;
unsigned char flag_huffman = 1;
unsigned char flag_wavelet_tree = 0;
unsigned char flag_entries = 1;
unsigned char flag_use_gpu = 1;

//global program parametrrs
unsigned int genome_length;
 struct FMIndex_WT *FM_index_WT;
 struct wavelet_tree *WT_root;
unsigned char max_bits = sizeof(unsigned long long int)*8;
unsigned char*genome;
unsigned int*count_table;
unsigned int*entries;


void error_handler(char err[], int code) {
  if(code != CL_SUCCESS) {
    printf("%s, Error code:%d\n", err, code);
    exit(EXIT_FAILURE);
  }
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
 struct compressedFMIndex *compressed_FM_index = NULL;
 unsigned int result_map;

 FILE *fp;
 FILE *fh_patterns;
 char *pattern = (char*)malloc(sizeof(char)*MAX_READ_LENGTH);
 /*
 for (i=0;i<argc;i++)
 {
  printf("%s\n",argv[i]);
 }*/

 //genome = load_genome_from_file_by_chunks(CHUNK_SIZE,filename_text,&genome_length);
 genome = load_genome_from_file(filename_text,&genome_length);

 printf("loaded genome\n");

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
  //free(genome); 
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
  
  char c;
  unsigned int *result = (unsigned int*) malloc (sizeof(unsigned int)*2);
  if (result == NULL)
  {
   printf("error pri alokacii\n");
   exit(-1);
  }

  printf("...approximate searching...\n");
  k = 0;
  count = 0;

  clock_t begin1 = clock();
  //get name of the read
  while (fgets(pattern, MAX_READ_LENGTH, fh_patterns) != NULL) 
  { 
   for (i=0;i<3;i++)
   {
    fgets(&pattern[i*READS_CHUNK], READS_CHUNK+2, fh_patterns );
   }

   pattern[strlen(pattern)-1]='\0';
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
  

  free(FM_index_WT->sampleSA);
  free(count_table);
  free(WT_root);
  free(FM_index_WT);
 }
 else if (flag_entries)
 {

  printf("rebuilding index...\n");
  fflush(stdout);

  rebuild_FM_index_into_entries(suffix_array,bwt);
  
  if (flag_use_gpu)
  {

    char *patterns_batch = (char *) malloc (sizeof(char) * (WARP_COUNT * NUM_OF_READS_PER_WARP * 64 + 1));
    i = 0;
    while (i != (WARP_COUNT*NUM_OF_READS_PER_WARP))
    {
      fgets(&patterns_batch[i*64], 70, fh_patterns); //header of pattern
      i++;
    }

    printf("%s\n",patterns_batch);
/* define platform */
  cl_platform_id platformID;
  ret = clGetPlatformIDs(1, &platformID, NULL);
  error_handler("Get platform error", ret); 

  cl_device_id deviceID;
  ret = clGetDeviceIDs(platformID, CL_DEVICE_TYPE_GPU, 1, &deviceID, NULL);
  error_handler("Get deviceID error", ret);

  cl_char vendorName[1024] = {0};
  cl_char deviceName[1024] = {0};
  ret = clGetDeviceInfo(deviceID, CL_DEVICE_VENDOR, sizeof(vendorName), vendorName, NULL);
  ret = clGetDeviceInfo(deviceID, CL_DEVICE_NAME, sizeof(deviceName), deviceName, NULL);
  printf("Connecting to %s %s...\n", vendorName, deviceName);

  /* define context */
  cl_context context;
  context = clCreateContext(NULL, 1, &deviceID, NULL, NULL, &ret);
  error_handler("define context error", ret);

  /* define command queue */
  cl_command_queue commandQueue;
  commandQueue = clCreateCommandQueue(context, deviceID, 0, &ret);
  error_handler("Create command queue error", ret);


  /* create memory objects */
  cl_mem inputMemObj;
  inputMemObj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(char)*WARP_COUNT*NUM_OF_READS_PER_WARP*64, NULL, &ret);
  error_handler("Create input buffer failed", ret);
  ret = clEnqueueWriteBuffer(commandQueue, inputMemObj, CL_TRUE, 0, sizeof(char)*WARP_COUNT*NUM_OF_READS_PER_WARP*64, (const void*)patterns_batch, 0, NULL, NULL);
  error_handler("Write to input buffer failed", ret);

  /* create memory objects */
  cl_mem inputMemObj_count_table;
  inputMemObj_count_table = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned int)*(alphabet_size+2), NULL, &ret);
  error_handler("Create input buffer failed", ret);
  ret = clEnqueueWriteBuffer(commandQueue, inputMemObj_count_table, CL_TRUE, 0, sizeof(unsigned int)*(alphabet_size+2), (const void*)count_table, 0, NULL, NULL);
  error_handler("Write to input buffer failed", ret);

  cl_mem inputMemObj_fm_index;
  inputMemObj_fm_index = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned int)*((genome_length+255) / 256 * 32), NULL, &ret);
  error_handler("Create input buffer failed", ret);
  ret = clEnqueueWriteBuffer(commandQueue, inputMemObj_fm_index, CL_TRUE, 0, sizeof(unsigned int)*((genome_length+255) / 256 * 32), (void*)entries, 0, NULL, NULL);
  error_handler("Write to input buffer failed", ret);

  // output need only be an array of inputSize / DEVICE_GROUP_SIZE long
  cl_mem outputMemObj;
  outputMemObj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(unsigned int)*WARP_COUNT*NUM_OF_READS_PER_WARP, NULL, &ret);
  error_handler("Create output buffer failed", ret);


  /* read kernel file from source */
  char fileName[] = "./sum.cl";
  char *sourceStr;  
  size_t sourceSize;

  FILE *file_source = fopen(fileName, "r");
  if(!file_source) {
    puts("Failed to load kernel file");
    exit(EXIT_FAILURE);
  }
  sourceStr = (char*)malloc(MAX_SOURCE_SIZE);
  sourceSize = fread(sourceStr, 1, MAX_SOURCE_SIZE, file_source);
  fclose(file_source);


  /* create program object */
  cl_program program = clCreateProgramWithSource(context, 1, (const char**)&sourceStr, (const size_t*)&sourceSize, &ret); 
  error_handler("create program failure", ret);
  ret = clBuildProgram(program, 1, &deviceID, NULL, NULL, NULL);
  if(ret != CL_SUCCESS) {
    puts("Build program error");
    size_t len;
    char buildLog[2048];
    clGetProgramBuildInfo(program, deviceID, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, &len);
    printf("%s\n", buildLog);
  }

  /* create kernel */
  cl_kernel kernel = clCreateKernel(program, "sum", &ret);
  error_handler("Create kernel failure", ret);

  unsigned int PATTERN_LENGTH = 64;
  /* set kernel arguments */
  ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&inputMemObj);
  error_handler("Set arg 1 failure", ret);
  ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&outputMemObj);
  error_handler("Set arg 2 failure", ret);
  ret = clSetKernelArg(kernel, 2, sizeof(int), (void *)&PATTERN_LENGTH);
  error_handler("Set arg 3 failure", ret);
  ret = clSetKernelArg(kernel, 3, sizeof(unsigned int) * 2 * NUM_OF_READS_PER_WARP, NULL);
  error_handler("Set arg 4 failure", ret);
  ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&inputMemObj_count_table);
  error_handler("Set arg 5 failure", ret);
  ret = clSetKernelArg(kernel, 5, sizeof(unsigned int) * 5, NULL);
  error_handler("Set arg 6 failure", ret);
  ret = clSetKernelArg(kernel, 6, sizeof(char) * NUM_OF_READS_PER_WARP * PATTERN_LENGTH, NULL);
  error_handler("Set arg 7 failure", ret);
  ret = clSetKernelArg(kernel, 7, sizeof(unsigned int) * 4 * NUM_OF_READS_PER_WARP, NULL); //indexes
  error_handler("Set arg 8 failure", ret);
  ret = clSetKernelArg(kernel, 8, sizeof(unsigned int) * DEVICE_GROUP_SIZE, NULL); //fm_index entry
  error_handler("Set arg 9 failure", ret);
  ret = clSetKernelArg(kernel, 9, sizeof(unsigned int) * DEVICE_GROUP_SIZE * 2, NULL); //count_table_results
  error_handler("Set arg 10 failure", ret);
  ret = clSetKernelArg(kernel, 10, sizeof(unsigned int) * DEVICE_GROUP_SIZE, NULL); //bitcounts
  error_handler("Set arg 11 failure", ret);
  ret = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *)&inputMemObj_fm_index);
  error_handler("Set arg 12 failure", ret);


   clock_t begin1 = clock();

  /* enqueue and execute */
  const size_t globalWorkSize = 768;
  const size_t localWorkSize = DEVICE_GROUP_SIZE;
  printf("global %lu, local %lu\n",globalWorkSize, localWorkSize);
  ret = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL);
  error_handler("Enqueue/execute failure", ret);

  unsigned int * results = (unsigned int*) malloc(sizeof(unsigned int)*NUM_OF_READS_PER_WARP*WARP_COUNT);

  ret = clEnqueueReadBuffer(commandQueue, outputMemObj, CL_TRUE, 0, sizeof(unsigned int)*NUM_OF_READS_PER_WARP*WARP_COUNT, results, 0, NULL, NULL);
  error_handler("Read output buffer fail", ret);


  clock_t end1 = clock();
  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("It took %lf seconds\n", time_spent1);

  }

  else
  {

  //char c;
  unsigned int *result = (unsigned int*) malloc (sizeof(unsigned int)*2);
  if (result == NULL)
  {
   printf("error pri alokacii\n");
   exit(-1);
  }

  printf("...approximate searching...\n");
  k = 0;
  count = 0;

  clock_t begin1 = clock();
  //get name of the read
  /*while (fgets(pattern, MAX_READ_LENGTH, fh_patterns) != NULL) 
  { 
   for (i=0;i<3;i++)
   {
    fgets(&pattern[i*READS_CHUNK], READS_CHUNK+2, fh_patterns );
   }
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
  }*/
  result_map = approximate_search_in_FM_index_entry("TCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCA",result);

  clock_t end1 = clock();
  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  printf("It took %lf seconds\n", time_spent1);
  printf("Total aligned reads: %d\n",count);
  printf("Total reads: %d\n",k);
  }
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
