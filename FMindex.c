#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FMindex.h"
#include "wavelet_tree.h"
#include <time.h>



//create FM Index which uses wavelet tree to compress BWT string
struct FMIndex_WT*build_FM_index_WT(unsigned int *suffix_array, unsigned char *bwt)
{
 unsigned int i;
 struct FMIndex_WT *FM_index_WT = (struct FMIndex_WT*) malloc(sizeof(struct FMIndex_WT));
 printf("...constructing huffman shaped wavelet tree... \n");
 //create count table and also create frequency table for huffman tree
 count_table = create_count_table(bwt);
 
 unsigned int*frequencies = (unsigned int*)malloc(sizeof(unsigned int)*strlen(alphabet));
 clock_t begin = clock();
 
 /*for (i=0;i<strlen(alphabet)-1;i++)
  frequencies[i] = FM_index_WT->count_table[i+1] - FM_index_WT->count_table[i];
 frequencies[i] = genome_length - FM_index_WT->count_table[i];*/

 for (i=0;i<strlen(alphabet);i++)
 {
  frequencies[i] = 1;
 }
 //store bwt in Huffman-shaped wavelet tree
 WT_root = build_huffman_shaped_WT(bwt,frequencies);
  clock_t end = clock();
 double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
 printf("construction of wavelet tree took %lf seconds\n", time_spent);
 fflush(stdout);

 printf("suffix array sample memory: %lu bytes\n",((genome_length+sample_SA_size)/sample_SA_size)*sizeof(unsigned int));
 FM_index_WT->end = find_end(suffix_array);
 FM_index_WT->sampleSA = create_sample_SA(suffix_array);
 free(suffix_array);
 printf("construction of FM Index finished\n");
 fflush(stdout);
 return FM_index_WT;
}

//create base FM Index
struct FMIndex*build_FM_index(unsigned int *suffix_array, unsigned char *bwt)
{
 unsigned char bits_per_char = get_min_bits_per_char(alphabet);
 struct FMIndex *FM_index = (struct FMIndex*) malloc(sizeof(struct FMIndex));
 FM_index->bwt = bwt;
 FM_index->alphabet = alphabet;
 FM_index->end = find_end(suffix_array);
 FM_index->sampleSA = create_sample_SA(suffix_array);
 free(suffix_array);
 FM_index->count_table = create_count_table(bwt);
 FM_index->occurence_table = create_occurence_table(bwt,alphabet);
 FM_index->alphabetically_encoded  = 0;
 alphabet_encode(FM_index->bwt,alphabet,genome_length);
 FM_index->alphabetically_encoded = 1;

 return FM_index;
}

//create FM Index using compressed BWT by MTF,RLE,Huffman
struct compressedFMIndex*build_compressed_FM_index(unsigned int *suffix_array, unsigned char *bwt, 
  unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman, unsigned int block_size)
{
 int i, count=0;
 struct compressedFMIndex *compressed_FM_index = (struct compressedFMIndex*) malloc(sizeof(struct compressedFMIndex));
 compressed_FM_index->block_size = block_size;
 compressed_FM_index->alphabet = alphabet;
 compressed_FM_index->flag_mtf = flag_mtf;
 compressed_FM_index->flag_runs = flag_runs;
 compressed_FM_index->flag_huffman = flag_huffman;
 compressed_FM_index->length = genome_length;
 compressed_FM_index->end = find_end(suffix_array);
 compressed_FM_index->sampleSA = create_sample_SA(suffix_array);
 free(suffix_array);
 compressed_FM_index->count_table = create_count_table(bwt);

 clock_t begin = clock();
 compressed_FM_index->array_of_blocks = compress_FMIndex(block_size,flag_mtf,flag_runs,flag_huffman,alphabet,bwt,&compressed_FM_index->length);
 clock_t end = clock();
 double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
 printf("trvanie kompresie: %lf", time_spent);
 printf("pocet blkov je %d\n",compressed_FM_index->length);
 fflush(stdout);
 
 //huffman_decode(compressed_FM_index->array_of_blocks[0].bitvector,compressed_FM_index->array_of_blocks[0].huffman_tree,compressed_FM_index->array_of_blocks[0].block_size);
 /*count_occ_in_compressed_FMIndex(38, 'G',compressed_FM_index);
 count_occ_in_compressed_FMIndex(39, 'G',compressed_FM_index);
 count_occ_in_compressed_FMIndex(40, 'G',compressed_FM_index);
 count_occ_in_compressed_FMIndex(41, 'G',compressed_FM_index);*/
 //printf("vysledok je %d",get_SA_value_compressed(1, genome_length, sample_SA_size, alphabet, 
  //compressed_FM_index->sampleSA,compressed_FM_index->count_table,compressed_FM_index->block_size, compressed_FM_index->array_of_blocks, flag_mtf, flag_runs, flag_huffman));
 //search_pattern_in_compressed_FM_index(compressed_FM_index, "GTTA",flag_mtf, flag_runs, flag_huffman);
 
 /*clock_t begin2 = clock();
 unsigned char *dcp = decompress_FMIndex(compressed_FM_index);
 clock_t end2 = clock();
 double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
 printf("trvanie DEkompresie: %lf", time_spent2);

 printf("DEcompression succesful\n");
 */
 
 for (i=0;i<=compressed_FM_index->length;i++){
  count = count + (compressed_FM_index->array_of_blocks[i].bitvector_length+7)/8;
  printf("bitvector length of %d block is %d\n",i,compressed_FM_index->array_of_blocks[i].bitvector_length);
 }
 printf("bwt is stored in %d bytes, and has %d characters\n",count,genome_length);
 return compressed_FM_index;
}

//find "end" of BWT string, it means position of zero rotation in BWT
unsigned int find_end(unsigned int *suffix_array)
{
unsigned int i = 0;
while(suffix_array[i]!=0)
 i++;
return i;
}

//procedure to sample input Suffix array according to input sample size
unsigned int *create_sample_SA(unsigned int *suffix_array)
{
 int i = 0;
 int sampled_index=0;
 int samples_count = (genome_length+sample_SA_size)/sample_SA_size;
 unsigned int *sampleSA = (unsigned int *) malloc((samples_count)*sizeof(unsigned int)); 
 if (sampleSA==NULL)
  {
   printf("Error. No memory when allocating sample of suffix array\n");
   exit(1);
  }
 while (i!=genome_length)
 {
  if (i%sample_SA_size == 0)
   sampleSA[sampled_index++]=suffix_array[i];
  i++;
 }
 return sampleSA; 
}

//procedure to get from sample suffix array  SA value corresponding to BWT position
//when input position not in sample suffix array, use LFmapping until returned position in sample suffix array
unsigned int get_SA_value(unsigned int bwt_position, unsigned char c, struct FMIndex *fm_index)
{
 unsigned int count = 0;
 unsigned char character;
 if (fm_index->alphabetically_encoded)
  character = c;
 else
  character = get_alphabet_index(fm_index->alphabet,c);
 while(bwt_position%sample_SA_size!=0)
 {
  c = fm_index->bwt[bwt_position];
  bwt_position = last_to_first_encoded(c, bwt_position, fm_index->alphabet,fm_index->count_table,fm_index->bwt,fm_index->occurence_table);
  count++;
 }
 unsigned int result = fm_index->sampleSA[bwt_position/sample_SA_size]+count;
 if (result>genome_length)
  result = result - genome_length;
 return result;
}

//get original position from position in BWT
unsigned int get_SA_value_compressed(unsigned int bwt_position, unsigned int length, unsigned char *alphabet, unsigned int*sampleSA, 
  unsigned int *count_table,unsigned int block_size, struct compressed_block*block, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman)
{
 unsigned int count = 0;
 unsigned int index = bwt_position/block_size;
 unsigned int position;
 unsigned char c;

 //until position of stored SA value is found, keep LFmapping position
 while(bwt_position%sample_SA_size!=0)
 {
  position = bwt_position - index*block_size;
  unsigned char *bwt = decompress_block(block[index].bitvector_length,block[index].bitvector,flag_mtf,flag_runs,flag_huffman,position+1,alphabet,block[index].huffman_tree);
  c = bwt[position];
  free(bwt);
  bwt_position = last_to_first_in_compressed_FMIndex(c,bwt_position,alphabet,count_table,block, block_size, flag_mtf,flag_runs,flag_huffman);
  index = bwt_position/block_size;
  count++;
 }
 unsigned int result = sampleSA[bwt_position/sample_SA_size]+count;
 if (result>=length)
  result = result - length;
 return result;
}

//procedure to get from sample suffix array  SA value corresponding to BWT position
//when input position not in sample suffix array, use LFmapping until returned position in sample suffix array
unsigned int get_SA_value_WT(unsigned int bwt_position)
{
 unsigned int count = 0;
 unsigned char c = wt_access(bwt_position, WT_root);
 while(bwt_position%sample_SA_size!=0)
 {
  c = wt_access(bwt_position, WT_root);
  bwt_position = last_to_first_WT(c,bwt_position,WT_root,count_table);
  count++;
 }
 unsigned int result = FM_index_WT->sampleSA[bwt_position/sample_SA_size]+count;
 if (result>=genome_length)
  result = result - genome_length;
 return result;
}

//procedure of LF-mapping
unsigned int last_to_first(unsigned char c, unsigned int bwt_position, unsigned char *alphabet,unsigned int *count_table,unsigned char*bwt,unsigned int**occurence_table)
{
 unsigned char character = get_alphabet_index(alphabet,c);
 unsigned int last = count_table[character] + count_occ(bwt,occurence_table,bwt_position,c,character);
 return last;
}

//procedure of LF-mapping in compressed BWT
unsigned int last_to_first_in_compressed_FMIndex(unsigned char c, unsigned int bwt_position, unsigned char*alphabet, unsigned int* count_table,struct compressed_block *block, unsigned int block_size, 
  unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman)
{
 unsigned int character = get_alphabet_index(alphabet,c);
 unsigned int last = count_table[character] + count_occ_in_compressed_FMIndex(block,block_size,bwt_position,c,flag_mtf, flag_runs, flag_huffman, alphabet);
 return last;
}

//procedure of LF-mapping in BWT string encoded to indexes in alphabet
unsigned int last_to_first_encoded(unsigned char c, unsigned int bwt_position, unsigned char *alphabet,unsigned int *count_table,unsigned char*bwt,unsigned int**occurence_table)
{
 unsigned int last = count_table[c] + count_occ(bwt,occurence_table,bwt_position,c,c);
 return last;
}

//procedure of LF-mapping in BWT string as WT
unsigned int last_to_first_WT(unsigned char c, unsigned int bwt_position, struct wavelet_tree *wtree_root,unsigned int *count_table)
{
 unsigned int character = get_alphabet_index(alphabet,c);
 unsigned int last = count_table[character] + wt_rank(c, bwt_position, wtree_root);
 return last;
}

//procedure which creates count table of FM Index for input string and alphabet
unsigned int *create_count_table(unsigned char *s)
{
 unsigned int alphabet_size = strlen(alphabet);
 unsigned int *count_table = (unsigned int *) calloc (alphabet_size+1,sizeof(unsigned int));
 int i = 0, j = 0, k = 0;
 for (i=0;i<genome_length;i++)
 {
  for (j=0;j<alphabet_size;j++)
  {
   if (s[i]==alphabet[j])
   {
    count_table[j]++;
    break;
   }
  }
 }
 j = count_table[0];
 count_table[0] = 0;
 for (i=1;i<alphabet_size;i++)
 {
  k = count_table[i];
  count_table[i] = count_table[i-1]+j;
  j = k;
 }

 count_table[alphabet_size] = genome_length;
 return count_table;
}

//procedure which creates sampled 2D occurence table for input string and alphabet
unsigned int **create_occurence_table(unsigned char *s, unsigned char *alphabet)
{
 unsigned int alphabet_size = strlen(alphabet);
 unsigned int samples_count = (genome_length+sample_OCC_size)/sample_OCC_size;
 unsigned int i,j;
 unsigned int index = 0;
 unsigned int **occurence_table = (unsigned int **)malloc(alphabet_size*sizeof(unsigned int*));
 for (i=0;i<alphabet_size;i++)
 {
  occurence_table[i] = (unsigned int *)malloc((samples_count+1)*sizeof(unsigned int));
  if (occurence_table[i]==NULL)
  {
   printf("Error. No memory when allocating occurence table\n");
   exit(1);
  }
 }

 //initialize number of occurences at zero position
 for (i=0;i<alphabet_size;i++)
 {
   occurence_table[i][0] = 0;
 }

 for (i=0;i<genome_length;i++)
 {
  index = (i+sample_OCC_size-1)/sample_OCC_size;
  if (i%sample_OCC_size==1)
   for (j=0;j<alphabet_size;j++)
    occurence_table[j][index]=occurence_table[j][index-1];
  for (j=0;j<alphabet_size;j++)
   if (s[i]==alphabet[j])
    occurence_table[j][index]++;     
 } 
 return occurence_table;
}

//procedure of FM Index for getting count of occurences of char c upto input position
//procedure uses sampled occurence table
unsigned int count_occ(unsigned char *s, unsigned int **occurence_table, unsigned int position, unsigned char c, unsigned char character)
{
 unsigned int count = 0;
 unsigned int a = position/sample_OCC_size;
 unsigned int bound = position - (a * sample_OCC_size);
 while (bound>0)
 {
  bound--;
  position--;
  if (s[position]==c)
   count++;
 }
 if (s[position]==c)
   count--;
return occurence_table[character][position/sample_OCC_size]+count;
}

//procedure of FM Index for getting count of occurences of char c upto input position when BWT is compressed in blocks
unsigned int count_occ_in_block(struct compressed_block *block, unsigned int position,unsigned char c,unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman, unsigned char*alphabet)
{
 unsigned int i;
 unsigned int count = 0;
 unsigned int alphabet_index = get_alphabet_index(alphabet,c);
 unsigned int occurences = block->occurences[alphabet_index];
 unsigned char *bwt = decompress_block(block->bitvector_length,block->bitvector,flag_mtf,flag_runs,flag_huffman,position,alphabet,block->huffman_tree);
 for (i=0;i<position;i++)
  if (bwt[i] == c)
    count++;
 free(bwt);
 return (occurences+count);
}

unsigned int count_occ_in_compressed_FMIndex(struct compressed_block *block, unsigned int block_size, unsigned int position, unsigned char c, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman, unsigned char*alphabet)
{
 unsigned int index = position/block_size;
 unsigned int remainder = position - index*block_size; //faster than modulo ?
 return count_occ_in_block(&block[index],remainder,c,flag_mtf, flag_runs, flag_huffman, alphabet);
}

unsigned int count_occ_in_decompressed_FMIndex(struct compressed_block *block, unsigned char*bwt,unsigned int block_size, unsigned int position, unsigned char c, unsigned char*alphabet)
{
  unsigned int index = position/block_size;
  unsigned int start = index*block_size;
  int count = 0;
  while (start<position)
  {
    if (bwt[start]==c)
      count++;
    start++;
  }
  if (bwt[position]==c)
      count++;
  unsigned char a = get_alphabet_index(alphabet,c);
  return (block[index].occurences[a]+count);
}

//procedur which inplace decompresses BWT of FM Index
unsigned char *decompress_FMIndex(struct compressedFMIndex *compressed_FM_index)
{
 unsigned int i, index = 0;
 unsigned char *decompressed = NULL;
 unsigned char*alphabet = compressed_FM_index->alphabet;
 unsigned char flag_mtf = compressed_FM_index->flag_mtf;
 unsigned char flag_runs = compressed_FM_index->flag_runs;
 unsigned char flag_huffman = compressed_FM_index->flag_huffman;
 unsigned int block_size = compressed_FM_index->block_size;
 struct compressed_block *cb = NULL;
 unsigned char *bwt = (unsigned char*) malloc (1);
 for (i=0;i<=compressed_FM_index->length;i++)
 {
  cb = &compressed_FM_index->array_of_blocks[i];
  decompressed = decompress_block(cb->bitvector_length, cb->bitvector, flag_mtf, flag_runs, flag_huffman, cb->block_size, alphabet, cb->huffman_tree);
  bwt = (unsigned char *) realloc(bwt,index+cb->block_size+1);
  strncpy(&bwt[index],decompressed,cb->block_size);
  free(decompressed);
  index = index + cb->block_size;
 free(compressed_FM_index->array_of_blocks[i].bitvector); 
 }
 bwt[index]='\0';
 return 0;
}

//procedure used for reversing BWT which is compressed
/*unsigned char* reverseBWT_compressed(unsigned char*bwt, unsigned int end, unsigned int* count_table, unsigned char*alphabet, 
  unsigned int block_size, struct compressed_block *block, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman)
{
  unsigned wow = end;
 unsigned int i=0;
 unsigned int j = genome_length-1;
 unsigned char a;
 unsigned char character;
 unsigned char *reversed = (unsigned char *)malloc(genome_length);
 if (reversed==NULL)
 {
  printf("Error. No memory when allocating memory when reversing string\n");
  exit(1);
 }
 while(i++!=genome_length)
 {
  a = bwt[end];
  reversed[j] = a;
  j--;
  //printf("end=%d,a je %d=%c\n",end,a,a);
  character = get_alphabet_index(alphabet,a);
  //printf("pricom %d + %d\n",count_table[character],count_occ_in_decompressed_FMIndex(block,bwt,block_size,end,a,alphabet));
  end = count_table[character]-1 + count_occ_in_decompressed_FMIndex(block,bwt,block_size,end,a,alphabet);
  if (end>genome_length){
    //printf("end je %d, length je %d\n",end,length);
    end = end - genome_length;
  }
  //printf("nova pozicia je %d\n",end);
 }
 reversed[genome_length]='\0';
 //printf("REVERSED JE\n%s\n",reversed);
 //free(bwt);
 return reversed;
}*/

//procedure which decode string of numbers to string of chars according to indexes of alphabet
void alphabet_decode (unsigned char *s, unsigned char * alphabet)
{
 unsigned int string_length = strlen(s);
 unsigned int i;
 for (i=0;i<string_length;i++)
 {
  s[i] = alphabet[s[i]];
 }
}

//procedure used for reversing BWT
unsigned char *reverseBWT(unsigned int end, unsigned char*alphabet,unsigned int*count_table,unsigned char*bwt, unsigned int**occurence_table, unsigned char alphabetically_encoded)
{
 int i = 0;
 int j = genome_length -1;
 char a;
 char *reversed = (char *)malloc(sizeof(char)*genome_length);
 if (reversed==NULL)
 {
  printf("Error. No memory when allocating memory when reversing string\n");
  exit(1);
 }
 if (alphabetically_encoded == 0)
 {
  while (i++!=genome_length)
  {
   a = bwt[end];
   reversed[j] = a;
   j--;
   end = last_to_first(a,end,alphabet,count_table,bwt,occurence_table);
  }
 }
 else
 {
  while (i++!=genome_length)
  {
   a = bwt[end];
   reversed[j] = a;
   j--;
   end = last_to_first_encoded(a,end,alphabet,count_table,bwt,occurence_table);
  }
  alphabet_decode(reversed,alphabet);
 }
return reversed;
}

//aux. procedure for getting index of char in alphabet
unsigned char get_alphabet_index(char *alphabet, unsigned char c)
{
 return (unsigned int)(strchr(alphabet,c) - alphabet);
}

void print_occurence_table(struct FMIndex *fm_index)
{
 int i,j;
 int alphabet_size = strlen(fm_index->alphabet);
 printf("Printing occurence table:\n");
 for (j=0;j<(genome_length+sample_OCC_size)/sample_OCC_size;j++)
 {
  printf("%d: ",j*sample_OCC_size);
  for (i=0;i<alphabet_size;i++)
  {
  printf("%d ",fm_index->occurence_table[i][j]);
  }
  printf("\n");
 }
}

void print_sample_SA(struct FMIndex *fm_index)
{
 int i;
 int samples = (genome_length+sample_SA_size)/sample_SA_size;
 int *sample_SA = fm_index->sampleSA;
 printf("Printing sample suffix array:\n");
 for (i = 0; i<samples;i++)
  printf("%d = %d\n",i*sample_SA_size,sample_SA[i]);
}

void print_count_table(struct FMIndex *fm_index)
{
 unsigned int i = 0;
 printf("Printing count table: \n");
 for (i = 0; i<strlen(fm_index->alphabet); i++)
  printf("%4c ",fm_index->alphabet[i]);
 printf("\n");
 for (i = 0; i<strlen(fm_index->alphabet); i++)
  printf("%4d ",fm_index->count_table[i]);
 printf("\n");
}

void print_info_fm_index(struct FMIndex *fm_index)
{
 printf("Length of BWT is: %d\n",genome_length);
 printf("Stored BWT is: %s--\n",fm_index->bwt);
 print_bit_vector(fm_index->bwt, fm_index->bitvector_length);
 printf("Reversed   is: %s--\n",reverseBWT(fm_index->end,fm_index->alphabet,fm_index->count_table,fm_index->bwt, fm_index->occurence_table,fm_index->alphabetically_encoded));
 print_count_table(fm_index);
 print_sample_SA(fm_index);
 print_occurence_table(fm_index);
}

//procedure for searching pattern in BWT of FM Index
unsigned int*search_pattern(struct FMIndex *fm_index, char *pattern)
{
 unsigned int*result=(unsigned int*)malloc(2*sizeof(unsigned int));
 unsigned int *count_table = fm_index->count_table;
 unsigned int pattern_length = strlen(pattern);
 int j = pattern_length - 1;
 int alphabet_index = get_alphabet_index(alphabet,pattern[j]);
 result[0] = count_table[alphabet_index];
 if (alphabet_index!=alphabet_size)
  result[1] = count_table[alphabet_index + 1]-1;
 else 
  result[1] = genome_length-1;
 while (j>0 && result[0]<result[1])
 {
  j--;
  alphabet_index = get_alphabet_index(alphabet,pattern[j]);
  result[0] = count_table[alphabet_index] + count_occ(fm_index->bwt,fm_index->occurence_table,result[0],alphabet_index,alphabet_index);
  result[1] = count_table[alphabet_index] + count_occ(fm_index->bwt,fm_index->occurence_table,result[1],alphabet_index,alphabet_index);
  }
 if (j>0 && result[0]>=result[1])
  result[1]--;
 return result;
}

//procedure for searching pattern in compressed BWT of FM Index
unsigned int*search_pattern_in_compressed_FM_index(struct compressedFMIndex *compressed_fm_index, char *pattern,unsigned char flag_mtf,unsigned char flag_runs, unsigned char flag_huffman)
{
 unsigned int*result=(unsigned int*)malloc(2*sizeof(unsigned int));
 unsigned int *count_table = compressed_fm_index->count_table;
 unsigned int pattern_length = strlen(pattern);
 unsigned int block_size = compressed_fm_index->block_size;
 int j = pattern_length - 1;
 int alphabet_index = get_alphabet_index(alphabet,pattern[j]);
 result[0] = count_table[alphabet_index];
 if (alphabet_index!=alphabet_size)
  result[1] = count_table[alphabet_index + 1]-1;
 else 
 {
  result[1] = genome_length-1;
 }
 while (j>0 && result[0]<result[1])
 {
  j--;
  alphabet_index = get_alphabet_index(alphabet,pattern[j]);
  
  result[0] = count_table[alphabet_index] + count_occ_in_compressed_FMIndex(compressed_fm_index->array_of_blocks,block_size,result[0],pattern[j],flag_mtf,flag_runs,flag_huffman,alphabet);
  result[1] = count_table[alphabet_index] + count_occ_in_compressed_FMIndex(compressed_fm_index->array_of_blocks,block_size,result[1],pattern[j],flag_mtf,flag_runs,flag_huffman,alphabet);
 }
 if (j>0 && result[0]>=result[1])
  result[1]--;
 return result;
}

//procedure for searching pattern in compressed BWT of FM Index
unsigned char search_pattern_in_FM_index_WT(char *pattern, unsigned char current_pattern_length, unsigned int*result)
{
 int j = current_pattern_length - 1;
 int alphabet_index = get_alphabet_index(alphabet,pattern[j]);
 unsigned int start,end;
 result[0] = count_table[alphabet_index];
 if (alphabet_index!=alphabet_size)
  result[1] = count_table[alphabet_index + 1]-1;
 else 
 {
  result[1] = genome_length-1;
 }
 while (j>0 && result[0]<result[1])
  {
   alphabet_index = get_alphabet_index(alphabet,pattern[--j]);
   result[0] = count_table[alphabet_index] + wt_rank(pattern[j], result[0], WT_root);
   result[1] = count_table[alphabet_index] + wt_rank(pattern[j], result[1], WT_root);
  }
 if (result[0]>=result[1])
 {
  result[1]--;
  return j;
 }
 return j;
}

//procedure which breaks input pattern and search each separately in FM Index
unsigned int*approximate_search(struct FMIndex *fm_index, unsigned char *pattern)
{
 unsigned int pattern_length = strlen(pattern);
 unsigned int div_pattern_length = (pattern_length+1)/(max_error+1);
 unsigned int j;
 int i; 
 unsigned int *result;
 unsigned char *patterns[max_error+1];
 for (i=0;i<max_error;i++)
 {
  patterns[i] = (unsigned char *)malloc(sizeof(unsigned char)*(div_pattern_length+1)); 
  memcpy( patterns[i], &pattern[i*div_pattern_length], div_pattern_length);
  patterns[i][div_pattern_length]= '\0';
 } 
 j = pattern_length-i*div_pattern_length;
 patterns[i] = (unsigned char *)malloc(sizeof(unsigned char)*(j+1)); 
 memcpy( patterns[i], &pattern[i*div_pattern_length], j);
 patterns[i][j]= '\0';

 for (i=max_error;i>=0;i--)
 {
  result = search_pattern(fm_index,patterns[i]);
 }
return result;
}

//procedure which breaks input pattern and search each separately in compressed FM index
unsigned int*approximate_search_in_compressed_FM_index(struct compressedFMIndex *compressed_fm_index, unsigned char *pattern, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman)
{
 unsigned int pattern_length = strlen(pattern);
 unsigned int div_pattern_length = (pattern_length+1)/(max_error+1);
 unsigned int j;
 int i; 
 unsigned int *result;
 unsigned char *patterns[max_error+1];
 for (i=0;i<max_error;i++)
 {
  patterns[i] = (unsigned char *)malloc(sizeof(unsigned char)*(div_pattern_length+1)); 
  memcpy( patterns[i], &pattern[i*div_pattern_length], div_pattern_length);
  patterns[i][div_pattern_length]= '\0';
 } 
 j = pattern_length-i*div_pattern_length;
 patterns[i] = (unsigned char *)malloc(sizeof(unsigned char)*(j+1)); 
 memcpy( patterns[i], &pattern[i*div_pattern_length], j);
 patterns[i][j]= '\0';

 for (i=max_error;i>=0;i--)
 {
  result = search_pattern_in_compressed_FM_index(compressed_fm_index,patterns[i],flag_mtf,flag_runs,flag_huffman);
 }
return result;
}

//procedure which breaks input pattern and search each separately in FM index - Wavelet tree
long long int approximate_search_in_FM_index_WT(unsigned char *pattern, unsigned int*result)
{
 unsigned int pattern_length = strlen(pattern);
 unsigned int k,temp, length,j;
 unsigned char pattern_half[pattern_length+2];
 unsigned int result_length, current_pattern_start;
 char current_error = 0;
 
 //search last part of the pattern
 unsigned char searching_part[pattern_length];
 unsigned int div_pattern_length = (pattern_length+1)/(max_error+1);
 unsigned int current_pattern_length = pattern_length;

 int i = max_error;
 current_pattern_start = pattern_length;
 while(i>=0 && current_error<=max_error)
 {
  if (div_pattern_length<current_pattern_start)
  {
  current_pattern_start = current_pattern_start - div_pattern_length;
  current_pattern_length = div_pattern_length;
  }
  else
  {
    current_pattern_length = current_pattern_start;
  current_pattern_start = 0;
  }

  memcpy( searching_part, &pattern[current_pattern_start], current_pattern_length);
  searching_part[current_pattern_length] = '\0';
  result_length = search_pattern_in_FM_index_WT(searching_part,current_pattern_length, result);
  if (!(result_length))
  {
   if (result[1]-result[0]<THRESHOLD)
   {
     
     while (result[0]!=result[1])
     {
      temp = get_SA_value_WT(result[0]) - current_pattern_start;
      result[0]++;
      /*unsigned char pattern_test[pattern_length+1];
      memcpy(pattern_test, &pattern[current_pattern_start], pattern_length-current_pattern_length);
      pattern_test[pattern_length-current_pattern_length] = '\0';*/
      k = 0;
      while (genome[temp]!=pattern[0])
      {
       temp++;
       k++;
      }
      if (k<=max_error+1)
      {
       unsigned char genome_test[pattern_length+2];
       if (k>0)
       {
        memcpy(genome_test, &genome[temp], pattern_length+1);
        genome_test[pattern_length+1] = '\0';
       }
       else 
       {
        memcpy(genome_test, &genome[temp], pattern_length+1);
        genome_test[pattern_length+1] = '\0';
       }
       if (align(pattern,genome_test))
        return temp;
      }
     }
    }
   }

   i--;
  }
return 0;
}

long long int approximate_search_gpu(unsigned char *patterns, unsigned int *results, unsigned int pattern_length, unsigned int reads_count)
{
  unsigned int k,temp, i;
  unsigned int current_pattern_start = 32;
  unsigned char pattern_test[100];
  unsigned int counter = 0;
  char o = 0;

  for (i = 0; i<reads_count*4;i++,i++)
  {
    o = 0;
    if (results[i+1]-results[i]<THRESHOLD)
    {
      while (results[i]!=results[i+1])
      {
        temp = get_SA_value_entry(results[i]);
        results[i]++;
        if (temp>=current_pattern_start)
        {
          temp = temp - current_pattern_start;
       
          memcpy(pattern_test, &patterns[i/4*pattern_length], pattern_length);
          pattern_test[pattern_length] = '\0';
          k = 0;
          while (genome[temp]!=pattern_test[0] && (k<=max_error+1))
          {
           temp++;
           k++;
          }
          if (k<=max_error+1)
          {
           unsigned char genome_test[pattern_length+2];
           if (k>0)
           {
            memcpy(genome_test, &genome[temp], pattern_length+1);
            genome_test[pattern_length+1] = '\0';
           }
           else 
           {
            memcpy(genome_test, &genome[temp], pattern_length+1);
            genome_test[pattern_length+1] = '\0';
          }
           if (align(pattern_test,genome_test))
           {
            counter++;
            printf("%d-%s-\n",i,pattern_test);
            printf("Approximate position: %u\n",temp);
            o = 1;
            break;
           }
          }
        }
      }
    }

    current_pattern_start = 0;
    i++;i++;
    if (o==0 && results[i+1]-results[i]<THRESHOLD)
    {
      while (results[i]!=results[i+1])
      {
        temp = get_SA_value_entry(results[i]);
        results[i]++;
        if (temp>=current_pattern_start)
        {
          temp = temp - current_pattern_start;
       
          memcpy(pattern_test, &patterns[i/4*pattern_length], pattern_length);
          pattern_test[pattern_length] = '\0';
          k = 0;
          while (genome[temp]!=pattern_test[0] && (k<=max_error+1))
          {
           temp++;
           k++;
          }
          if (k<=max_error+1)
          {
           unsigned char genome_test[pattern_length+2];
           if (k>0)
           {
            memcpy(genome_test, &genome[temp], pattern_length+1);
            genome_test[pattern_length+1] = '\0';
           }
           else 
           {
            memcpy(genome_test, &genome[temp], pattern_length+1);
            genome_test[pattern_length+1] = '\0';
          }
           if (align(pattern_test,genome_test))
           {
            counter++;
            printf("%d-%s-\n",i,pattern_test);
            printf("Approximate position: %u\n",temp);
            break;
           }
          }
        }
      }
    }


  }
return counter;
  

}
//procedure which breaks input pattern and search each separately in FM index - Wavelet tree
long long int approximate_search_in_FM_index_entry(unsigned char *pattern, unsigned int*result)
{
 unsigned int pattern_length = strlen(pattern);
 unsigned int k,temp, length,j;
 unsigned char pattern_half[pattern_length+2];
 unsigned int result_length, current_pattern_start;
 char current_error = 0;
 
 //search last part of the pattern
 unsigned char searching_part[pattern_length];
 unsigned int div_pattern_length = (pattern_length+1)/(max_error+1);
 unsigned int current_pattern_length = pattern_length;

 int i = max_error;
 current_pattern_start = pattern_length;

 while(i>=0)
 {
  
  if (div_pattern_length<current_pattern_start)
  {
  current_pattern_start = current_pattern_start - div_pattern_length;
  current_pattern_length = div_pattern_length;
  }
  else
  {
    current_pattern_length = current_pattern_start;
  current_pattern_start = 0;
  }

  memcpy( searching_part, &pattern[current_pattern_start], current_pattern_length);
  searching_part[current_pattern_length] = '\0';
  result_length = search_pattern_in_FM_index_entry(searching_part,current_pattern_length, result);
  if (!(result_length))
  {
   if (result[1]-result[0]<THRESHOLD)
   {
    while (result[0]!=result[1])
    {

      temp = get_SA_value_entry(result[0]);
      result[0]++;
      if (temp>=current_pattern_start)
      {
       temp = temp - current_pattern_start;
       
       /*unsigned char pattern_test[pattern_length+1];
       memcpy(pattern_test, &pattern[current_pattern_start], pattern_length-current_pattern_length);
       pattern_test[pattern_length-current_pattern_length] = '\0';*/
       k = 0;
       while (genome[temp]!=pattern[0] && (k<=max_error+1))
       {
        temp++;
        k++;
       }
       if (k<=max_error+1)
       {
        unsigned char genome_test[pattern_length+2];
        if (k>0)
        {
         memcpy(genome_test, &genome[temp], pattern_length+1);
         genome_test[pattern_length+1] = '\0';
        }
        else 
        {
         memcpy(genome_test, &genome[temp], pattern_length+1);
         genome_test[pattern_length+1] = '\0';
        }
        if (align(pattern,genome_test))
         return temp;
       }
     }
    }
   }
  }
  i--;
 }
return 0;
}

//procedure for calculating square of dynamic programming matrix
unsigned int calculate(int i, int j, int k)
{
 unsigned int ret_value;
 if (i>j)
  if (i>k)
   ret_value = i;
  else 
   ret_value = k;
 else
  if (j>k) 
   ret_value = j;
  else
   ret_value = k;
 if (ret_value<0)
  return 0;
 else
  return ret_value;
}

//procedure for getting score of 2 characters of 2 strings
unsigned int score(char a, char b)
{
  if (a==b)
    return 1;
  else 
    return 0;
}

//procedure which returns max value in array
unsigned int get_max_array(unsigned char* array, unsigned int length)
{
 unsigned int i, max = 0, maxIndex = 0;
 for (i=0;i<length;i++) 
 {
  if (array[i]>max)
  {
    max = array[i];
    maxIndex = i;
  }
 }
 return maxIndex;
}

//procedure which returns direction of previous value in dynamic programming matrix
unsigned char get_nearest_max(unsigned char i, unsigned char j, unsigned char k)
{
 if (i>=j)
  if (i>=k)
   return 1;
  else 
   return 3;
 else
  if (j>=k) 
   return 2;
  else
   return 3;
}

//procedure for aligning two strings up to max error
unsigned char align(char *p1, char*p2)
{
 
 unsigned int corrects;
 int i,j;
 unsigned int p1length = strlen(p1), p2length = strlen(p2);
 unsigned char matrix[p1length+1][p2length+1];

 for (i=0;i<p1length+1;i++)
  for (j=0;j<p2length+1;j++)
    matrix[i][j] = 0;

 for (i=0;i<=max_error;i++)
  matrix[0][i] = 0;
 for (i=0;i<=max_error;i++)
  matrix[i][0] = 0;
 
 for (i=1;i<=max_error;i++){
  for (j=1;j<=i+max_error;j++){
   matrix[i][j] = calculate(matrix[i-1][j-1] + score(p1[i-1],p2[j-1]),matrix[i-1][j] - 1,matrix[i][j-1] - 1);
  }
}

 for (i=1+max_error;i<=p1length;i++){
  for (j=i-max_error;j<=i+max_error;j++){
    matrix[i][j] = calculate(matrix[i-1][j-1] + score(p1[i-1],p2[j-1]),matrix[i-1][j] - 1,matrix[i][j-1] - 1);
   }
}

 corrects = get_max_array(&matrix[p1length][p1length-max_error],2*max_error+1);
 corrects = corrects+p1length-max_error;

 for (i=p1length-max_error;i<=p1length;i++)
 {
  if ((i-matrix[i][p2length])<=max_error){
    return 1;
  }
 }
  
  for (i=p2length-max_error;i<=p2length;i++)
  {
    if ((i-matrix[p1length][i])<=max_error)
    {
      return 1;
    }
  }
  return 0;
}

void rebuild_FM_index_into_entries(unsigned int*suffix_array, unsigned char*bwt)
{
 unsigned int i = 0, a;
 unsigned int chars_in_entry = 256;
  
  count_table = create_count_table(bwt);
  unsigned int bitword; unsigned int bit_index;

  unsigned int sample_sa_part = chars_in_entry/sample_SA_size;  //8*4 = 32
  unsigned int bits_part = chars_in_entry/8/sizeof(unsigned int); //8*4 = 32
  unsigned int children_part = bits_part + 1; //9*4 = 36
  unsigned int aux_part = 2; //2*4 = 8
  unsigned int occ_root_part = 2; //2*4 = 8
  unsigned int total = (genome_length+255) / 256;
  
  entries = (unsigned int*) calloc((total+1)*32,sizeof(unsigned int)); //128/4 = 32

  printf("For entries was totally allocated %u bytes\n", (total)*32*4);
  //iterators
  unsigned int sa_iter = 0;
  unsigned int bit_iter = 0;
  unsigned int left_iter = 0;
  unsigned int right_iter = 0;
  unsigned int root_counter_iter = 0;
  unsigned int left_counter_iter = 0;
  unsigned int right_counter_iter = 0;


  unsigned int left_count = 0;
  unsigned int right_count = 0;
  unsigned int temp_count;
  unsigned int WARP_SIZE = 128;
  unsigned int remainder;

 unsigned int frequencies[4];
 for (i=0;i<strlen(alphabet);i++)
  frequencies[i] = 1;

 //store bwt in Huffman-shaped wavelet tree
 struct wavelet_tree *root = build_huffman_shaped_WT(bwt,frequencies);
 struct wavelet_tree *left = root->left;
 struct wavelet_tree *right = root->right; 
 printf("suffix array sample memory: %lu bytes\n",((genome_length+sample_SA_size)/sample_SA_size)*sizeof(unsigned int));
 unsigned int*sampleSA = create_sample_SA(suffix_array);
  
 i = 0;
 while (--total)
 {

  //store wt root - positions 8..15
  //printf("zacinam zapisovat %d ROOT na poziciu %d\n",bit_iter,i); //8 pre pozicie 0..7
 right_count = 0;
 for (a = 0; a<bits_part/2;a++)
 {
  entries[i++] = root->bitvector[bit_iter]>>32; 
  entries[i++] = root->bitvector[bit_iter];
  right_count += __builtin_popcountll(root->bitvector[bit_iter]);
  bit_iter++;
 }
 left_count = 256 - right_count;
 
 //printf("zacinam zapisovat LEFT_INFO na poziciu %d\n",i); //1 pre info? ..8
 //store children info about how many
 entries[i++] = (left_count+32)/32; 
 
 
 //printf("zacinam zapisovat %d left CHILDREN %d na poziciu %d\n",left_iter, left_count, i); // 9..17
 //store wt_children positions
 bitword = 0;
 bit_index = 0;
 for (a = (i/32)*256; a<((i/32)+1)*256;a++)
 {
  if (bwt[a]=='A' || bwt[a]=='T')
  {
   if (bwt[a]=='T')
   {
    bitword = bitword<<1;
    bitword = bitword + 1;
    bit_index++;
   }
   else
   {
    bitword = bitword<<1;
    bit_index++;
   }

   if (bit_index==32)
   {
    //printf("ukadam bitword %u\n",bitword);
    entries[i++] = bitword;
    bit_index = 0;
    bitword = 0;
   }
  }
 }

if (bit_index){
 entries[i++] = bitword << (32-bit_index);
 //printf("ukadama bitword\n");
}
else
  entries[i++] = 0;

//printf("zacinam zapisovat %d rigthCHILDREN %d na poziciu %d\n",right_iter, right_count, i);
 bitword = 0;
 bit_index = 0;
 for (a = (i/32)*256; a<((i/32)+1)*256;a++)
 {
  if (bwt[a]=='C' || bwt[a]=='G')
  {
   if (bwt[a]=='C')
   {
    bitword = bitword<<1;
    bitword = bitword + 1;
    bit_index++;
   }
   else
   {
    bitword = bitword<<1;
    bit_index++;
   }

   if (bit_index==32)
   {
    //printf("ukadam bitword %u\n",bitword);
    entries[i++] = bitword;
    bit_index = 0;
    bitword = 0;
   }
  }
 }

if (bit_index) {
 entries[i++] = bitword << (32-bit_index);
 //printf("ukadama bitword\n");
 }
 //printf("zacinam zapisovat OCC_ROOT na poziciu %d\n",i);  //2hodnoty 18..19
 //store counters
 
 temp_count = 0;
  for (a = 0; a < 2; a++)
  {
    temp_count += __builtin_popcount (entries[i+a-18]);
  }
  entries[i++] = temp_count;//count_unset_bits(root->bitvector,root_counter_iter-1,root->bitcount_table);
  for (a = 2; a < 5; a++)
  {
    temp_count += __builtin_popcount (entries[i+a-19]);
  }
  entries[i++] = temp_count;//count_set_bits(root->bitvector,root_counter_iter-1,root->bitcount_table);


 /*printf("zacinam zapisovat %d left_OCC_CHILDREN na poziciu %d\n",left_counter_iter,i); //2hodnoty 20..21
 fflush(stdout);*/
 if (left_counter_iter)
 {
  entries[i++] = count_unset_bits(left->bitvector,left_counter_iter-1,left->bitcount_table);
  entries[i++] = count_set_bits(left->bitvector,left_counter_iter-1,left->bitcount_table);
 }
 else
 {
  entries[i++] = 0;
  entries[i++] = 0;
 }
 left_counter_iter = left_counter_iter + left_count;
 
 //printf("zacinam zapisovat %d right_OCC_CHILDREN na poziciu %d\n",right_counter_iter,i); //2hodnoty 22..23
 if (right_counter_iter)
 {
 entries[i++] = count_unset_bits(right->bitvector,right_counter_iter-1,right->bitcount_table);
 entries[i++] = count_set_bits(right->bitvector,right_counter_iter-1,right->bitcount_table);
 }
 else
 {
  entries[i++] = 0;
  entries[i++] = 0;
 }
 right_counter_iter = right_counter_iter + right_count;

 /*printf("zacinam zapisovat %d SA_VALS na poziciu %d\n",sa_iter,i); //8hodnot 24..31
 fflush(stdout);*/
 //storing SA_values
 for (a = 0; a<sample_sa_part;a++)
  entries[i++] = sampleSA[sa_iter++];

}

//printf("idem zapisat zvysky\n");
//printf("zacinam zapisovat %d ROOT na poziciu %d\n",bit_iter,i); //8 pre pozicie 0..7
 right_count = 0;
 for (a = 0; a<bits_part/2;a++)
 {
  if (bit_iter*64 >= genome_length)
  {
    i++;i++;
  }
  else
  {
   entries[i++] = root->bitvector[bit_iter]>>32; 
   entries[i++] = root->bitvector[bit_iter];
   right_count += __builtin_popcountll(root->bitvector[bit_iter]);
   bit_iter++;
  }
 }
 left_count = (genome_length - (i/32)*256) - right_count;
 
 //printf("zacinam zapisovat LEFT_INFO %d na poziciu %d\n",((left_count+32)/32),i); //1 pre info? ..8
 //store children info about how many
 entries[i++] = (left_count+32)/32; 
 
 //printf("zacinam zapisovat %d left CHILDREN %d na poziciu %d\n",left_iter, left_count, i); // 9..17
 //store wt_children positions
 bitword = 0;
 bit_index = 0;
 for (a = (i/32)*256; a<genome_length;a++)
 {
  if (bwt[a]=='A' || bwt[a]=='T')
  {
   if (bwt[a]=='T')
   {
    bitword = bitword<<1;
    bitword = bitword + 1;
    bit_index++;
   }
   else
   {
    bitword = bitword<<1;
    bit_index++;
   }

   if (bit_index==32)
   {
    //printf("ukadam bitword %u\n",bitword);
    entries[i++] = bitword;
    bit_index = 0;
    bitword = 0;
   }
  }
 }
if (bit_index){
 entries[i++] = bitword << (32-bit_index);
 //printf("ukadama bitword\n");
}
else
  entries[i++] = 0;

//printf("zacinam zapisovat %d rigthCHILDREN %d na poziciu %d\n",right_iter, right_count, i);
 bitword = 0;
 bit_index = 0;
 for (a = (i/32)*256; a<genome_length;a++)
 {
  if (bwt[a]=='C' || bwt[a]=='G')
  {
   if (bwt[a]=='C')
   {
    bitword = bitword<<1;
    bitword = bitword + 1;
    bit_index++;
   }
   else
   {
    bitword = bitword<<1;
    bit_index++;
   }

   if (bit_index==32)
   {
    //printf("ukadam bitword %u\n",bitword);
    entries[i++] = bitword;
    bit_index = 0;
    bitword = 0;
   }
  }
 }

if (bit_index) {
 entries[i++] = bitword << (32-bit_index);
 //printf("ukadama bitword\n");
 }

 while (i%32 != 18)
  i++;

//store counters

temp_count = 0;
  for (a = 0; a < (genome_length - (i/32)*256)/32; a++)
  {
    temp_count += __builtin_popcount (entries[i+a-18]);
  }
  entries[i++] = temp_count;//count_unset_bits(root->bitvector,root_counter_iter-1,root->bitcount_table);
  for (a = 2; a < (genome_length - (i/32)*256)/32; a++)
  {
    temp_count += __builtin_popcount (entries[i+a-19]);
  }
  entries[i++] = temp_count;//count_set_bits(root->bitvector,root_counter_iter-1,root->bitcount_table);

/*fflush(stdout);
 exit(-1);*/

 /*printf("zacinam zapisovat %d left_OCC_CHILDREN na poziciu %d\n",left_counter_iter,i); //2hodnoty 20..21
 fflush(stdout);*/
 if (left_counter_iter)
 {
  entries[i++] = count_unset_bits(left->bitvector,left_counter_iter-1,left->bitcount_table);
  entries[i++] = count_set_bits(left->bitvector,left_counter_iter-1,left->bitcount_table);
 }
 else
 {
  entries[i++] = 0;
  entries[i++] = 0;
 }
 left_counter_iter = left_counter_iter + left_count;
 
 //printf("zacinam zapisovat %d right_OCC_CHILDREN na poziciu %d\n",right_counter_iter,i); //2hodnoty 22..23
 if (right_counter_iter)
 {
 entries[i++] = count_unset_bits(right->bitvector,right_counter_iter-1,right->bitcount_table);
 entries[i++] = count_set_bits(right->bitvector,right_counter_iter-1,right->bitcount_table);
 }
 else
 {
  entries[i++] = 0;
  entries[i++] = 0;
 }
 right_counter_iter = right_counter_iter + right_count;

 /*printf("zacinam zapisovat %d SA_VALS na poziciu %d\n",sa_iter,i); //8hodnot 24..31
 fflush(stdout);*/
 //storing SA_values
 for (a = 0; a<sample_sa_part;a++)
  entries[i++] = sampleSA[sa_iter++];

free(sampleSA);
}


//procedure for searching pattern in compressed BWT of FM Index
unsigned char search_pattern_in_FM_index_entry(char *pattern, unsigned char current_pattern_length, unsigned int*result)
{

 int j = current_pattern_length - 1;
 unsigned char alphabet_index = get_alphabet_index(alphabet,pattern[j]);
 unsigned int start,end;
 
 result[0] = count_table[alphabet_index];
 result[1] = count_table[alphabet_index + 1]-1;

//printf("--------------\n");
 while (j>0 && result[0]<result[1])
  {
   alphabet_index = get_alphabet_index(alphabet,pattern[--j]);
   
   //printf("res: %d %d\n",result[0], result[1]);
   unsigned int entry_index1 = (result[0]>>8) * 32; //in each entry up to 256 chars
   unsigned int in_entry_index1 = result[0]&255; //v danom entry

   unsigned int entry_index2 = (result[1]>>8) * 32; //in each entry up to 256 chars
   unsigned int in_entry_index2 = result[1]&255; //v danom entry

  __builtin_prefetch (&entries[entry_index1], 0, 1);
  __builtin_prefetch (&entries[entry_index2], 0, 1);

   result[0] = count_table[alphabet_index] + wt_rank_entry(alphabet_index, result[0], entry_index1, in_entry_index1);
   result[1] = count_table[alphabet_index] + wt_rank_entry(alphabet_index, result[1], entry_index2, in_entry_index2);
  }


//printf("res: %d %d\n",result[0], result[1]);
 if (result[0]>=result[1])
 {
  result[1]--;
  return j;
 }
  
 return j;
}


unsigned int wt_rank_entry(unsigned char c, unsigned int position, unsigned int entry_index, unsigned int in_entry_index)
{
 unsigned int i;
 unsigned int count = 0;
 unsigned int result = 0;
 unsigned int a = in_entry_index;
 unsigned int remainder = in_entry_index;

  //spocitaju sa 1ky od zaciatku entry po dany in_entry_index
  for (i=0;i<in_entry_index/32;i++)
  {
    count = count + __builtin_popcount(entries[entry_index+i]);
    remainder = remainder - 32;
  }
  if (remainder)
   count = count + __builtin_popcount(entries[entry_index+i]>>(32-remainder));

 if (c == 0 || c == 3)
 {
  //bolo nutne spocitat nuly, preto invertujem
  count = in_entry_index - count;

  //spocitam 1 do pozicie count
  remainder = count;
  for (i=0;i<count/32;i++)
    {
     result = result + __builtin_popcount(entries[entry_index+9+i]);
     remainder = remainder - 32;
    }
  
   if (remainder)
    result = result + __builtin_popcount(entries[entry_index+9+i]>>(32-remainder));

  if (c==0)
   return (count - result + entries[entry_index + 20]);
  
  else
    return (result + entries[entry_index + 21]);
 }

 else 
 {
  a = entries[entry_index + 8]; //kolko bolo v left?
  remainder = count;
  for (i=0;i<count/32;i++)
  {
   result = result + __builtin_popcount(entries[entry_index+9+a+i]);
   remainder = remainder - 32;
  }

  if (remainder)
   result = result + __builtin_popcount(entries[entry_index+9+a+i]>>(32-remainder));

 if (c==2)
   return (count - result + entries[entry_index + 22]);

  else
   return (result + entries[entry_index + 23]);
 }
}

unsigned int get_SA_value_entry(unsigned int bwt_position)
{
 unsigned int count = 0;
 unsigned int entry_index; //in each entry up to 256 chars
 unsigned short int in_entry_index; //v danom entry
 unsigned int character;
 unsigned char c;
 while(bwt_position%sample_SA_size!=0)
 {
  entry_index = (bwt_position/256) * 32;
  in_entry_index = bwt_position%256;

  __builtin_prefetch (&entries[entry_index], 0, 3);
  c = wt_access_entry(entry_index, in_entry_index);
  character = get_alphabet_index(alphabet,c);
  bwt_position = count_table[character] + wt_rank_entry(character, bwt_position, entry_index, in_entry_index);
  count++;
 }
 
 entry_index = (bwt_position/256) * 32;
 in_entry_index = bwt_position%256;
 unsigned int result = entries[entry_index + 24 + in_entry_index/32] + count;
 
 if (result>=genome_length)
  result = result - genome_length;

 return result;
}


//function which returns character on input position
unsigned char wt_access_entry(unsigned int entry_index, unsigned short int in_entry_index)
{
 unsigned int i;
 unsigned int count = 0;
 unsigned int a = in_entry_index;
 unsigned int remainder = in_entry_index;
 //printf("indexy su %d %d \n",entry_index,in_entry_index);
  //spocitaju sa 1ky od zaciatku entry po dany in_entry_index

  for (i=0;i<in_entry_index/32;i++)
  {
    count = count + __builtin_popcount(entries[entry_index+i]);
    //printf("count je %d\n",count);
    remainder = remainder - 32;
  }
  if (remainder){
    //printf("hm remmainder %d\n",remainder);
   count = count + __builtin_popcount(entries[entry_index+i]>>(32-remainder));
  }
  


  //printf("wut count je %d, i je%d\n",count,i);
  if (get_bit_4byte(entries[entry_index + i],remainder))
  {

    //printf("pocitam 1, entry_index %d\n",entry_index);
    //aky znak na count pozicii?
    a = entries[entry_index+8];
    
    //remainder = count%32;

    //printf("pre count %u, pristupujem na poziciu %u %u\n",count,count/sample_SA_size,remainder);
    //printf("entries[%d] je %u\n",entry_index + 9 + a + count/32, entries[entry_index + 9 + a + count/32]);
    if (get_bit_4byte(entries[entry_index + 9 + a + count/32], count%32))
      return 'C';
    else
      return 'G';
  }

  else
  {
    
    //pre nulu je potrebne spocitat nuly, preto invertujem
    count = in_entry_index - count;

    if (get_bit_4byte(entries[entry_index + 9 + count / 32], count % 32))
      return 'T';
    else
      return 'A';
  }
}

