#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FMindex.h"
#include "compression.h"
#include <time.h>

unsigned int THRESHOLD = 3;

//create FM Index which uses wavelet tree to compress BWT string
struct FMIndex_WT*build_FM_index_WT(unsigned int *suffix_array, unsigned char *bwt)
{
 unsigned int i;
 struct FMIndex_WT *FM_index_WT = (struct FMIndex_WT*) malloc(sizeof(struct FMIndex_WT));
 printf("...constructing huffman shaped wavelet tree... \n");
 //create count table and also create frequency table for huffman tree
 FM_index_WT->count_table = create_count_table(bwt);
 unsigned int*frequencies = (unsigned int*)malloc(sizeof(unsigned int)*strlen(alphabet));
 clock_t begin = clock();
 for (i=0;i<strlen(alphabet)-1;i++)
  frequencies[i] = FM_index_WT->count_table[i+1] - FM_index_WT->count_table[i];
 frequencies[i] = genome_length - FM_index_WT->count_table[i];

 //store bwt in Huffman-shaped wavelet tree
 FM_index_WT->WT_root = build_huffman_shaped_WT(bwt,frequencies);
  clock_t end = clock();
 double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
 printf("construction of wavelet tree took %lf seconds\n", time_spent);
 fflush(stdout);

 printf("suffix array sample memory: %d bytes\n",((genome_length+sample_SA_size)/sample_SA_size)*sizeof(unsigned int));
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
unsigned int get_SA_value_WT(unsigned int bwt_position, struct FMIndex_WT *fm_index_wt)
{
 unsigned int count = 0;
 unsigned char c = wt_access(bwt_position, fm_index_wt->WT_root);
 while(bwt_position%sample_SA_size!=0)
 {
  //printf("******%d %c %d\n",bwt_position,c,count);
  c = wt_access(bwt_position, fm_index_wt->WT_root);
  //printf("**** access pozicie %d je %c\n",bwt_position,c);
  bwt_position = last_to_first_WT(c,bwt_position,fm_index_wt->WT_root,fm_index_wt->count_table);
  count++;
 }
 //printf("---pos %d RESULT: %d %d \n",bwt_position,fm_index_wt->sampleSA[bwt_position/fm_index_wt->sample_SA_size],count);
 unsigned int result = fm_index_wt->sampleSA[bwt_position/sample_SA_size]+count;
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
 //printf("LF %c %d %d %d %d\n",c,character,count_table[character], wt_rank(c, bwt_position, wtree_root),last);
 return last;
}

//procedure which creates count table of FM Index for input string and alphabet
unsigned int *create_count_table(unsigned char *s)
{
 unsigned int alphabet_size = strlen(alphabet);
 unsigned int *count_table = (unsigned int *) calloc (alphabet_size,sizeof(unsigned int));
 int i = 0;
 int j = 0;
 int k = 0;
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
 return count_table;
}

//procedure which creates sampled 2D occurence table for input string and alphabet
unsigned int **create_occurence_table(unsigned char *s, unsigned char *alphabet)
{
 unsigned int alphabet_size = strlen(alphabet);
 unsigned int samples_count = (genome_length+sample_OCC_size)/sample_OCC_size;
 unsigned int i;
 unsigned int j;
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
  //printf("dekopmresujem blok cislo %d\n",i);
  decompressed = decompress_block(cb->bitvector_length, cb->bitvector, flag_mtf, flag_runs, flag_huffman, cb->block_size, alphabet, cb->huffman_tree);
  bwt = (unsigned char *) realloc(bwt,index+cb->block_size+1);
  strncpy(&bwt[index],decompressed,cb->block_size);
  free(decompressed);
  index = index + cb->block_size;
  //free(compressed_FM_index->array_of_blocks[i].occurences);
  free(compressed_FM_index->array_of_blocks[i].bitvector); 
 }
 bwt[index]='\0';
 //printf("velkost bwt je %d\n",index);
 //printf("%s\n",bwt);
//printf("Reversing .... \n");
//return reverseBWT_compressed(bwt,index,compressed_FM_index->end,compressed_FM_index->count_table, alphabet,block_size,compressed_FM_index->array_of_blocks,flag_mtf,flag_runs,flag_huffman);
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
unsigned int get_alphabet_index(unsigned char *alphabet, unsigned char c)
{
 int i = 0;
 while(c!=alphabet[i])
  i++;
 return i;
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
 //printf("prvy znak %c index %d\n",pattern[j],alphabet_index);
 result[0] = count_table[alphabet_index];
 if (alphabet_index!=alphabet_size)
  result[1] = count_table[alphabet_index + 1]-1;
 else 
  result[1] = genome_length-1;
 //printf("rozsah je %d az %d \n",result[0],result[1]);
 while (j>0 && result[0]<result[1])
 {
  j--;
  alphabet_index = get_alphabet_index(alphabet,pattern[j]);
  //printf("%d-ty znak %c index %d\n",j,pattern[j],alphabet_index);
  result[0] = count_table[alphabet_index] + count_occ(fm_index->bwt,fm_index->occurence_table,result[0],alphabet_index,alphabet_index);
  result[1] = count_table[alphabet_index] + count_occ(fm_index->bwt,fm_index->occurence_table,result[1],alphabet_index,alphabet_index);
 
 /*int test = result[0];
 printf("-----------------------------\n");
 while (test<result[1])
 {
 printf("sa value %d je %d\n",test,get_SA_value(test,fm_index->bwt[result[0]],fm_index));
 test++;
 }
 */

 }
 if (j>0 && result[0]>=result[1])
  result[1]--;
 printf("Range of positions of pattern are %d to %d \n",result[0],result[1]);
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
 //printf("prvy znak %c index %d\n",pattern[j],alphabet_index);
 result[0] = count_table[alphabet_index];
 if (alphabet_index!=alphabet_size)
  result[1] = count_table[alphabet_index + 1]-1;
 else 
 {
  result[1] = genome_length-1;
 }
 //printf("\n****rozsah je %d az %d \n",result[0],result[1]);
 while (j>0 && result[0]<result[1])
 {
  j--;
  alphabet_index = get_alphabet_index(alphabet,pattern[j]);
  
  result[0] = count_table[alphabet_index] + count_occ_in_compressed_FMIndex(compressed_fm_index->array_of_blocks,block_size,result[0],pattern[j],flag_mtf,flag_runs,flag_huffman,alphabet);
  result[1] = count_table[alphabet_index] + count_occ_in_compressed_FMIndex(compressed_fm_index->array_of_blocks,block_size,result[1],pattern[j],flag_mtf,flag_runs,flag_huffman,alphabet);
 //printf("****rozsah je %d az %d \n",result[0],result[1]);
 /*int test = result[0];
 printf("-----------------------------\n");
 while (test<result[1])
 {
 printf("sa value %d je %d\n",test,get_SA_value(test,fm_index->bwt[result[0]],fm_index));
 test++;
 }
 */

/*
****rozsah je 5031 az 9975 
****rozsah je 11196 az 12427 
****rozsah je 2839 az 3145 
****rozsah je 723 az 804 
****rozsah je 5210 az 5232 
****rozsah je 1332 az 1338

****rozsah je 14928 az 19999 
****rozsah je 13669 az 14927 
****rozsah je 18409 az 18720 
****rozsah je 4623 az 4710 
****rozsah je 1182 az 1204 
****rozsah je 10256 az 10262
*/

 }
 if (j>0 && result[0]>=result[1])
  result[1]--;
 
 printf("  rozsah je %d az %d \n",result[0],result[1]);
 /*unsigned int test = result[0];
 while (test<result[1])
 {
  printf("*************************************\n");
 printf("hladam test %d, gl %d\n",test,compressed_fm_index->genome_length);
 printf("sa value %d je %d\n",test,get_SA_value_compressed(test, compressed_fm_index->genome_length,  compressed_fm_index->sample_SA_size, alphabet,compressed_fm_index->sampleSA,count_table,block_size, compressed_fm_index->array_of_blocks, flag_mtf, flag_runs, flag_huffman));
 test++;
 }*/
 return result;
}

//procedure for searching pattern in compressed BWT of FM Index
unsigned int*search_pattern_in_FM_index_WT(struct FMIndex_WT *FM_index_WT, char *pattern)
{
 struct wavelet_tree *wtree_root = FM_index_WT->WT_root;
 unsigned int*result = (unsigned int*)malloc(2*sizeof(unsigned int));
 unsigned int *count_table = FM_index_WT->count_table;
 unsigned int pattern_length = strlen(pattern);
 int j = pattern_length - 1;
 int alphabet_index = get_alphabet_index(alphabet,pattern[j]);
 //printf("prvy znak %c index %d\n",pattern[j],alphabet_index);
 result[0] = count_table[alphabet_index];
 if (alphabet_index!=alphabet_size)
  result[1] = count_table[alphabet_index + 1]-1;
 else 
 {
  result[1] = genome_length-1;
 }
 //printf("\n****rozsah je %d az %d \n",result[0],result[1]);
 while (j>0 && result[0]<=result[1])
 {
  --j;
  alphabet_index = get_alphabet_index(alphabet,pattern[j]);
  //printf("problem %d\n",alphabet_index);
  //fflush(stdout);
  result[0] = count_table[alphabet_index] + wt_rank(pattern[j], result[0], wtree_root);
  result[1] = count_table[alphabet_index] + wt_rank(pattern[j], result[1], wtree_root);
 
 //printf("**%d**rozsah je %d az %d \n",j,result[0],result[1]);
 /*int test = result[0];
 printf("-----------------------------\n");
 while (test<result[1])
 {
 printf("sa value %d je %d\n",test,get_SA_value(test,fm_index->bwt[result[0]],fm_index));
 test++;
 }
 */


 }
 if (j>0 && result[0]>=result[1]){
  printf("HEJ --\n");
  --result[1];
}
 
 //printf("  rozsah je %d az %d \n",result[0],result[1]);
 
 /*unsigned int test = result[0]-1;
 while (test<result[1])
 {
 //printf("*************************************\n");
 //printf("hladam test %d, gl %d\n",test,FM_index_WT->genome_length);
 printf("  sa value %d je %d\n",test,get_SA_value_WT(test, FM_index_WT));
 test++;
 }*/

 return result;
}

//procedure for searching pattern in compressed BWT of FM Index
unsigned int*threshold_search_pattern_in_FM_index_WT(struct FMIndex_WT *FM_index_WT, char *pattern, unsigned int *result_length, unsigned int *result)
{
 struct wavelet_tree *wtree_root = FM_index_WT->WT_root;
 //unsigned int*result = (unsigned int*)malloc(2*sizeof(unsigned int));
 unsigned int *count_table = FM_index_WT->count_table;
 unsigned int pattern_length = strlen(pattern);
 int j = pattern_length - 1;
 int alphabet_index = get_alphabet_index(alphabet,pattern[j]);
 //printf("prvy znak %c index %d\n",pattern[j],alphabet_index);
 result[0] = count_table[alphabet_index];
 if (alphabet_index!=alphabet_size)
  result[1] = count_table[alphabet_index + 1]-1;
 else 
 {
  result[1] = genome_length-1;
 }
 //printf("\n****rozsah je %d az %d \n",result[0],result[1]);
 while (j>0 && (result[1]-result[0])>THRESHOLD)
 {
  j--;
  alphabet_index = get_alphabet_index(alphabet,pattern[j]);
  
  result[0] = count_table[alphabet_index] + wt_rank(pattern[j], result[0], wtree_root);
  result[1] = count_table[alphabet_index] + wt_rank(pattern[j], result[1], wtree_root);
 //printf("**%d**rozsah je %d az %d \n",j,result[0],result[1]);
 /*int test = result[0];
 printf("-----------------------------\n");
 while (test<result[1])
 {
 printf("sa value %d je %d\n",test,get_SA_value(test,fm_index->bwt[result[0]],fm_index));
 test++;
 }
 */
 *result_length = j;

 }
 if (j>0 && result[0]>=result[1]){
  result[1]--;
}
 
 /*printf("  rozsah je %d az %d dlzka %d\n",result[0],result[1],pattern_length-j);
 
 unsigned int test = result[0];
 while (test<result[1])
 {
 //printf("*************************************\n");
 //printf("hladam test %d, gl %d\n",test,FM_index_WT->genome_length);
 printf("    sa value %d je %d\n",test,get_SA_value_WT(test, FM_index_WT));
 test++;
 }*/

 return result;
}

//procedure which breaks input pattern and search each separately in FM Index
unsigned int*approximate_search(struct FMIndex *fm_index, unsigned char *pattern)
{
 unsigned int pattern_length = strlen(pattern);
 unsigned int div_pattern_length = (pattern_length+1)/(max_error+1);
 unsigned int j;
 int i; 
 unsigned int *result;
 printf("%s dlzka vzoru je %d\n",pattern,pattern_length);
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
  printf("%d. vzor je %s\n",i,patterns[i]);
  result = search_pattern(fm_index,patterns[i]);
  printf("pozicie su %d az %d\n-------------------------------------\n",result[0],result[1]);
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
 printf("%s dlzka vzoru je %d\n",pattern,pattern_length);
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
  printf("%d. vzor je %s\n",i,patterns[i]);
  result = search_pattern_in_compressed_FM_index(compressed_fm_index,patterns[i],flag_mtf,flag_runs,flag_huffman);
  printf("pozicie su %d az %d\n",result[0],result[1]);
 }
return result;
}

//procedure which breaks input pattern and search each separately in FM index - Wavelet tree
long long int approximate_search_in_FM_index_WT(struct FMIndex_WT *FM_index_WT, unsigned char *pattern, unsigned int *result)
{
 unsigned int pattern_length = strlen(pattern);
 unsigned int div_pattern_length = (pattern_length+1)/(max_error+1);
 unsigned int result_length = 0;
 unsigned int j,k;
 int i; 
 //char *p1 = (unsigned char*)malloc(sizeof(unsigned char)*(pattern_length+max_error+1));
 unsigned int m,n,count,innerStart,innerEnd,results_count=0;
 int current_error = 0;
 //unsigned int *result;
 unsigned int *perfect_results = (unsigned int *)malloc(sizeof(unsigned int)*(max_error+1)*2*4);
 
 unsigned char pattern_half[pattern_length];
 unsigned int* SA_results = NULL;

 //printf("%s dlzka vzoru je %d\n",pattern,pattern_length);
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
 
 i = 0;
 k = 0;
 
 //druha verzia:
  //rozdelit na casti
  //kazdu cast rozdelit na polovice
    //kazdu polovicu vyhladat pomocou backward a forward
      //ak su vysledky rovnake = 0 chyb
      //ak sa vysledky lisia = 1 chyba
      //ak v jednej polovici nie su vysledky = 1 chyba
      //ak v obidvoch poloviciach nie su vysledky = 2 chyby
    //uchovat spravne vysledky: zaciatok, dlzka, vysledky
  // 

 while(i<=max_error && current_error<=max_error)
 {
  //search i-th part of pattern 
  //printf("%d. vzor je %s\n",i,patterns[i]);
  result = search_pattern_in_FM_index_WT(FM_index_WT,patterns[i]);

  //if search fails:
  //  1. increase number of errors
  //  2. break searching part into halves
  //  3. each half maximize upon threshold and store
  //else store result
  if (result[1]<=result[0]){
    current_error++;
    //search first half and 
    strncpy(pattern_half,patterns[i],strlen(patterns[i])/2);
    //printf("kopirujem %d znakov\n",strlen(patterns[i])/2);
    pattern_half[strlen(patterns[i])/2] = '\0';
    //printf("  error - prvy retazec je %s a dlzka %d\n",pattern_half,strlen(pattern_half));
    
    printf("    results: %d %d\n",result[0],result[1]);
    
    threshold_search_pattern_in_FM_index_WT(FM_index_WT,pattern_half,&result_length,result);
    
    printf("    results: %d %d\n",result[0],result[1]);
    //printf("dlzkya je %d\n",strlen(pattern_half)-result_length);
    
     perfect_results[k++]=i*div_pattern_length+result_length;
     perfect_results[k++]=strlen(pattern_half)-result_length;
     perfect_results[k++]=result[0];
     perfect_results[k++]=result[1];
    
    //printf("results count je %d, roz %d %d\n",results_count, result[1],result[0]);
  results_count = 1+results_count + (result[1]-result[0]);
    //printf("    results: %d %d\n",result[0],result[1]);
   // printf("  druhy retazec je %s\n",&patterns[i][strlen(patterns[i])/2]);
    threshold_search_pattern_in_FM_index_WT(FM_index_WT,&patterns[i][strlen(patterns[i])/2],&result_length,result);
    //printf("    results: %d %d\n",result[0],result[1]);
    //printf("dlzkya je %d\n",result_length);
     perfect_results[k++]=i*div_pattern_length+strlen(pattern_half)+result_length;
     perfect_results[k++]=strlen(patterns[i])-result_length-strlen(pattern_half);
     perfect_results[k++]=result[0];
     perfect_results[k++]=result[1];

     //printf("results count je %d, roz %d %d\n",results_count, result[1],result[0]);
  results_count = 1+results_count + (result[1]-result[0]);

  }
  else
  {
   perfect_results[k++]=i*div_pattern_length;
   perfect_results[k++]=strlen(patterns[i]);
   perfect_results[k++]=result[0];
   perfect_results[k++]=result[1];
  }
  //printf("results count je %d\n",results_count);
  results_count = 1+results_count + (result[1]-result[0]);
  //printf("pozicie su %d az %d, akt. err: %d\n-------------------------------------\n",result[0],result[1],current_error);
  i++;
 }
 
 j=0;
 //printf("***Combining phase ***\n");
 if (i>max_error && current_error<=max_error)
 {
  SA_results = (unsigned int *)malloc(sizeof(unsigned int)*(results_count));
  //printf("alokujem %d pre p.results a %d pre SA_res\n",((max_error+1)*2*3),results_count);
  //printf("errors: %d Results are:\n", current_error);
  for (i=0;i<k;i=i+4){
     //printf("%d. vysledok bez chyby s dlzkou %d je %d-%d: ",perfect_results[i],perfect_results[i+1],perfect_results[i+2],perfect_results[i+3]);
     //printf("%dcisla;",perfect_results[i+3]-perfect_results[i+2]);
     SA_results[j++] = perfect_results[i+3]-perfect_results[i+2];
     //printf("%d\n",perfect_results[i+3]-perfect_results[i+2]);
     while(perfect_results[i+2]!=perfect_results[i+3])
     {
      //count SA values
      //printf("%d,",get_SA_value_WT(perfect_results[i+2],FM_index_WT));
      SA_results[j++] = get_SA_value_WT(perfect_results[i+2],FM_index_WT)-perfect_results[i];
      //and store them
      perfect_results[i+2]++;
    }
    //printf(";\n");
  }
 }
 else {
  /*printf("READ NOT MAPPED - There were more than %d error(s).\n",max_error);
  fflush(stdout);*/
  free(SA_results);
 free(perfect_results);
 for (i=0;i<=max_error;i++)
  free(patterns[i]);
 return 0;

 }
 //printf("p.results is %d and sa_res is %d\n",k,j);
//fflush(stdout);
count = 0;
n = 0;
//printf("k je %d i je %d\n",k,i);
//k/4 je pocet iteracii
for (i=0;i<=max_error;i++){
 //nacitam pocet cisiel
 n = count;
 count = count + SA_results[count];
 //printf("%d cyklus ide od %d do %d\n",i,n,count);
 n++;
 while(n<=count)
 {
  //zvolim n-tu hodnotu
  current_error = i;
  //good_counter = 0;
  //printf("idem porovnavat %d v %d cykle od %d do %d\n",SA_results[n],i,n,count);
  innerStart = count+1;
  for(j=i+1;j<k/4;j++){ //nie je nutne pocitat pre posledne kedze pocet chyb
   innerEnd = innerStart+SA_results[innerStart];
   while (innerStart<innerEnd)
   {
    innerStart++;
    //printf("  porovnavam s %d od %d do %d\n",SA_results[innerStart],innerStart,innerEnd);

    if(SA_results[n]>SA_results[innerStart])
      m = SA_results[n]-SA_results[innerStart];
    else m = SA_results[innerStart]-SA_results[n];

    if ((m+current_error)<=max_error)
    {
     current_error += m-1;
     //good_counter++;
     //printf("zhoda s %d chybou\n",m);
     innerStart = innerEnd;
    }

   }
   if (++current_error>max_error){
    break;
  }
   innerStart++;

  }
  //printf("vyhladavanie skoncilo s %d chybaimi\n",current_error);
  if (current_error<=max_error){
   //printf("READ MAPPED - predpokladane riesenie je na urovni %d +- %d, pocet chyb je %d\n",SA_results[n],max_error,current_error);
   i = max_error+1;//endd main forloop
   
   if (SA_results!=NULL)
    free(SA_results);
   free(perfect_results);
   for (i=0;i<=max_error;i++)
   {
    free(patterns[i]);
   }
   
   return 1;
  }

  n++;
 }

 count++;
}
if (current_error>max_error)
 //printf("READ NOT MAPPED - no combination found\n");
/*
 for (i=0;i<j;i++){
  for (k=i+1;k<j;k++){
    if(SA_results[i]>SA_results[k])
      m = SA_results[i]-SA_results[k];
    else m = SA_results[k]-SA_results[i];

   printf("%d %d, m = %d\n",SA_results[i],SA_results[k],m);
   if (m<=max_error)
   {
    printf("predpokladane riesenie je na urovni %d +- %d\n",SA_results[i],max_error);
    //strncpy(p1,);
    //align(char *p1, char*p2, int error)

    return SA_results[i];
   }
  }
 }*/
 if (SA_results!=NULL)
  free(SA_results);
 free(perfect_results);
 for (i=0;i<=max_error;i++){
  free(patterns[i]);
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
    return -1;
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
void align(char *p1, char*p2, int error)
{
 unsigned char **matrix;
 unsigned int corrects;
 unsigned int alignment_length;
 char *alignment1;
 char *alignment2;
 int i,j;
 unsigned int nearestMax;
 unsigned int currentX, currentY;
 unsigned int p1length = strlen(p1), p2length = strlen(p2);
 
 //alignment = (char *) malloc((p1length+error+1)*sizeof(char));
 
 matrix = (unsigned char **) malloc(sizeof(unsigned char*) * (p1length + 1));
 for (i=0;i<=p1length;i++)
  matrix[i] = (unsigned char *) calloc(p2length + 1,sizeof(unsigned char)); 
 for (i=0;i<=error;i++)
  matrix[0][i] = 0;
 for (i=0;i<=error;i++)
  matrix[i][0] = 0;
 
 for (i=1;i<=error;i++){
  for (j=1;j<=i+error;j++){
   matrix[i][j] = calculate(matrix[i-1][j-1] + score(p1[i-1],p2[j-1]),matrix[i-1][j] - 1,matrix[i][j-1] - 1);
   printf(" %d ",matrix[i][j]);
  }
  printf("\n");
}
 for (i=1+error;i<=p1length;i++){
  for (j=i-error;j<=i+error;j++){
    matrix[i][j] = calculate(matrix[i-1][j-1] + score(p1[i-1],p2[j-1]),matrix[i-1][j] - 1,matrix[i][j-1] - 1);
    printf(" %d ",matrix[i][j]);
   }
  printf("\n");
}
 corrects = get_max_array(&matrix[p1length][p1length-error],2*error+1);
 corrects = corrects+p1length-error;
 if (p1length-error>matrix[p1length][corrects])
  printf("No alignment found\n");
 else
 {
   alignment_length = p1length+(p1length - matrix[p1length][corrects])+1;
   alignment1 = (char *) malloc(alignment_length*sizeof(char));
   alignment2 = (char *) malloc(alignment_length*sizeof(char));
   alignment1[alignment_length-1] = '\0';
   alignment2[alignment_length-1] = '\0';
   currentY = corrects;
   currentX = p1length;
   for (i = alignment_length-2;i>0;i--)
   {
    if (currentX==0 || currentY==0)
      break;
    nearestMax = get_nearest_max(matrix[currentX-1][currentY-1],matrix[currentX-1][currentY],matrix[currentX][currentY-1]);
    //if (matrix[currentX][currentY]-1 == matrix[currentX-1][currentY-1]) || (matrix[currentX][currentY]+1 == matrix[currentX-1][currentY-1])
    if (nearestMax==1){
      alignment1[i] = p1[currentX-1];
      alignment2[i] = p2[currentY-1];
      currentX = currentX - 1;
      currentY = currentY - 1; 
    }
    else if (nearestMax==2){
      //printf("nearestMax2 je %d\n",matrix[currentX-1][currentY]);
      if (matrix[currentX-1][currentY]-1 == matrix[currentX][currentY]){
       alignment1[i] = p1[currentX-1];
       alignment2[i] = '-';
       currentX = currentX - 1;
      }
      else
      {
       //printf("cant get there %d %d\n",currentX,currentY);
       alignment1[i] = p1[currentX-1];
       alignment2[i] = p2[currentY-1];
       currentX = currentX - 1;
       currentY = currentY - 1; 
      }
    }
    else{
      //printf("nearestMax3 je %d\n",matrix[currentX][currentY-1]);
      if (matrix[currentX][currentY-1]-1 == matrix[currentX][currentY]){
       alignment1[i] = '-';
       alignment2[i] = p2[currentY-1];
       currentY = currentY - 1;
      }
      else
      {
       //printf("cant get there %d %d\n",currentX,currentY);
       alignment1[i] = p1[currentX-1];
       alignment2[i] = p2[currentY-1];
       currentX = currentX - 1;
       currentY = currentY - 1; 
      }
    }

   }
   //printf("%d i je %d, c %d %d\n",nearestMax,i,currentX,currentY);
   if (currentX==1)
   {
    if (currentY==1){
     alignment1[i] = p1[0];
     alignment2[i] = p2[0];
    }
    else{
     alignment1[i] = p1[0];
     alignment2[i] = '-';
    }
   }
   else if (currentY==1)
   {
    alignment1[i] = '-';
    alignment2[i] = p2[0];
   } 
   while (i>=0){
    alignment1[i] = ' ';
    alignment2[i] = ' ';
    i--;
   }
   printf("%s\n",alignment1);
   printf("%s\n",alignment2);
   printf("errors: %d\n",alignment_length - p1length-1);
 }

 free(alignment1);
 free(alignment2);
 for (i=0;i<=p1length;i++)
  free(matrix[i]);
 free(matrix);
}
