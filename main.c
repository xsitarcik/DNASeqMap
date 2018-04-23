#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bwt.h"
#include "FMindex.h"
#include <stdint.h>
#include "compression.h"
#include <CL/cl.h>

#define DEVICE_GROUP_SIZE 32
#define NUM_OF_READS_PER_WARP 1
#define WARP_COUNT 192
#define MAX_SOURCE_SIZE 0x100000
//MAIN PARAMETERS:
//for constructing auxiliary tables of FM Index
unsigned int sample_OCC_size = 2; //in reality it's *64, MUST SET TO 2
unsigned int sample_SA_size = 32;
unsigned char max_error = 1;
unsigned int THRESHOLD = 500;

//>SRR493095.1 M00282:31:000000000-A0FFK:1:1:13945:1807 length=150



unsigned char load = 0;

//program parameters
unsigned char save = 0;
unsigned char *save_name;

unsigned char*filename_text;// = "alt_Celera_chr15.fa";
unsigned char *load_name;// = "alt_Celera_chr15_bwt_withoutN.txt";
unsigned char*filename_patterns;// = "SRR493095final.fasta";
//unsigned char*filename_patterns = "data/kratke_ready.fasta";

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
unsigned char flag_use_gpu = 0;

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
  printf(" Required input parameter - choose only one\n");
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
 struct compressedFMIndex *compressed_FM_index = NULL;
 unsigned int result_map;


 naive_compress("AAGCT");

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
  else if (strcmp(argv[i],"-c") == 0)
    sample_OCC_size = convert_string_to_int(argv[++i]);
  else if (strcmp(argv[i],"-w") == 0)
    flag_wavelet_tree = 1;
  else if (strcmp(argv[i],"-n") == 0)
    flag_entries = 1;
  else if (strcmp(argv[i],"-g") == 0)
    flag_use_gpu = 1;

  else if (strcmp(argv[i],"-r") == 0)
  {
    filename_text = (char *) realloc (filename_text, sizeof(char) * strlen(argv[++i])+1);
    strcpy(filename_text, argv[i]);
  }
  else if (strcmp(argv[i],"-p") == 0)
  {
    filename_patterns = (char *) realloc (filename_patterns,sizeof(char)*strlen(argv[++i])+1);
    strcpy(filename_patterns, argv[i]);
  }
  else if (strcmp(argv[i],"-i") == 0)
  {
    load = 1;
    load_name = (char *) realloc (load_name, sizeof(char)*strlen(argv[++i])+1);
    strcpy(load_name, argv[i]);
  }
  else if (strcmp(argv[i],"-b") == 0)
  {
    save = 1;
    save_name = (char *) realloc (save_name, sizeof(char)*strlen(argv[++i])+1);
    strcpy(save_name, argv[i]);
  }
}

 FILE *fp;
 FILE *fh_patterns;
 char *pattern = (char*)malloc(sizeof(char)*MAX_READ_LENGTH);
 
 if (!(load) && !(save))
 {
  genome = load_genome_from_file_by_chunks(CHUNK_SIZE,filename_text,&genome_length);
 //genome = load_genome_from_file(filename_text,&genome_length);
 printf("Successfully loaded genome with length %u\n",genome_length);
 }

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

 //build FMIndex using compression methods
 if (flag_compress)
 {
  compressed_FM_index = build_compressed_FM_index(suffix_array,bwt,flag_mtf, flag_runs, flag_huffman, block_size);

  /*int*result = approximate_search_in_compressed_FM_index(compressed_FM_index,"TAATCGGTGGGAGTATTCAACGTGATGAAGAC",flag_mtf,flag_runs,flag_huffman);
  printf("results: %d %d\n",result[0],result[1]);
  printf("ide sa na vypocet\n");
  fflush(stdout);*/
  
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

   printf("nacital: -%s-, %d\n",pattern,i);
   fflush(stdout);
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
  ret = clBuildProgram(program, 1, &deviceID, "-cl-mad-enable", NULL, NULL);
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
  ret = clSetKernelArg(kernel, 7, sizeof(unsigned int) * 8 * NUM_OF_READS_PER_WARP, NULL); //indexes
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
  const size_t globalWorkSize = 6144;
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
  clock_t begin1 = clock();
  //get name of the read
  while (1) 
  { 
   char line[256];
   pattern[0] = '\0';
   i = 0;
   while (1)
    {
     if (fgets(line, sizeof(line), fh_patterns)!= NULL);
      else{
        o = 1;
        break;
      }
     if (line[0] == '>')
      break;
     else {
      strcpy(&pattern[i],line);
      i = i + strlen(line) - 1;
    }

    }
   if (o)
    break;
   /*for (i=0;i<3;i++) CTCAGCTTAGGACCCGACTAACCCAGAGCGGACGAGCCTTCCTCTGGAAACCTTAGTCAATCGGTGGACG
   {
    fgets(&pattern[i*READS_CHUNK], READS_CHUNK+2, fh_patterns );
   }*/
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
