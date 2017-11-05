#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FMindex.h"
#include "compression.h"

struct FMIndex*build_FM_index(unsigned int *suffix_array, unsigned int sample_SA_size, unsigned int sample_OCC_size, unsigned int genome_length, unsigned char *bwt, unsigned char *alphabet)
{
  unsigned char bits_per_char = get_min_bits_per_char(alphabet);
  printf("%s\n",bwt);
 struct FMIndex *FM_index = (struct FMIndex*) malloc(sizeof(struct FMIndex));
 FM_index->bwt = bwt;
 FM_index->sample_SA_size = sample_SA_size;
 FM_index->sample_OCC_size = sample_OCC_size;
 FM_index->alphabet = alphabet;
 FM_index->length = genome_length;
 FM_index->end = find_end(suffix_array);
 FM_index->sampleSA = create_sample_SA(suffix_array,sample_SA_size,genome_length);
 free(suffix_array);
 FM_index->count_table = create_count_table(bwt,genome_length,alphabet);
 FM_index->occurence_table = create_occurence_table(bwt,genome_length,alphabet,sample_OCC_size);
 FM_index->alphabetically_encoded  = 0;
 alphabet_encode(FM_index->bwt,alphabet,genome_length);
 FM_index->alphabetically_encoded = 1;
 //FM_index->bwt = bit_pack(FM_index->bwt,genome_length,bits_per_char,&FM_index->bitvector_length);
 count_occ(FM_index->bwt,FM_index->occurence_table,41,2,2,sample_OCC_size);
 return FM_index;
}

struct compressedFMIndex*build_compressed_FM_index(unsigned int *suffix_array, unsigned int sample_SA_size, unsigned int sample_OCC_size, unsigned int genome_length, unsigned char *bwt, unsigned char *alphabet, 
  unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman, unsigned int block_size)
{
  printf("%s\n",bwt);
 struct compressedFMIndex *compressed_FM_index = (struct compressedFMIndex*) malloc(sizeof(struct compressedFMIndex));
 compressed_FM_index->sample_SA_size = sample_SA_size;
 compressed_FM_index->block_size = block_size;
 compressed_FM_index->alphabet = alphabet;
 compressed_FM_index->flag_mtf = flag_mtf;
 compressed_FM_index->flag_runs = flag_runs;
 compressed_FM_index->flag_huffman = flag_huffman;
 compressed_FM_index->length = genome_length;
 compressed_FM_index->genome_length = genome_length;
 compressed_FM_index->end = find_end(suffix_array);
 compressed_FM_index->sampleSA = create_sample_SA(suffix_array,sample_SA_size,genome_length);
 free(suffix_array);
 compressed_FM_index->count_table = create_count_table(bwt,genome_length,alphabet);
 compressed_FM_index->array_of_blocks = compress_FMIndex(block_size,flag_mtf,flag_runs,flag_huffman,alphabet,bwt,&compressed_FM_index->length);
 printf("pocet blkov je %d\n",compressed_FM_index->length);
 //huffman_decode(compressed_FM_index->array_of_blocks[0].bitvector,compressed_FM_index->array_of_blocks[0].huffman_tree,compressed_FM_index->array_of_blocks[0].block_size);
 /*count_occ_in_compressed_FMIndex(38, 'G',compressed_FM_index);
 count_occ_in_compressed_FMIndex(39, 'G',compressed_FM_index);
 count_occ_in_compressed_FMIndex(40, 'G',compressed_FM_index);
 count_occ_in_compressed_FMIndex(41, 'G',compressed_FM_index);*/
 //printf("vysledok je %d",get_SA_value_compressed(1, genome_length, sample_SA_size, alphabet, 
  //compressed_FM_index->sampleSA,compressed_FM_index->count_table,compressed_FM_index->block_size, compressed_FM_index->array_of_blocks, flag_mtf, flag_runs, flag_huffman));
 //search_pattern_in_compressed_FM_index(compressed_FM_index, "GTTA",flag_mtf, flag_runs, flag_huffman);
 unsigned char *dcp = decompress_FMIndex(compressed_FM_index);
 //printf("%s\n",decompress_block(compressed_FM_index->array_of_blocks[0].bitvector_length,compressed_FM_index->array_of_blocks[0].bitvector,flag_mtf,flag_runs,flag_huffman,50,alphabet,compressed_FM_index->array_of_blocks[0].huffman_tree));
 return compressed_FM_index;
}

unsigned int find_end(unsigned int *suffix_array)
{
unsigned int i = 0;
while(suffix_array[i]!=0)
 i++;
return i;
}

unsigned int *create_sample_SA(unsigned int *suffix_array,unsigned int sample_size, unsigned int array_size)
{
 int i = 0;
 int sampled_index=0;
 int samples_count = (array_size+sample_size)/sample_size;
 unsigned int *sampleSA = (unsigned int *) malloc((samples_count)*sizeof(unsigned int)); 
 if (sampleSA==NULL)
  {
   printf("Error. No memory when allocating sample of suffix array\n");
   exit(1);
  }
 while (i!=array_size)
 {
  if (i%sample_size == 0)
   sampleSA[sampled_index++]=suffix_array[i];
  i++;
 }
 return sampleSA; 
}

//procedure to get SA value from sample suffix array
unsigned int get_SA_value(unsigned int bwt_position, unsigned char c, struct FMIndex *fm_index)
{
 //printf("%c znak\n",c);
 unsigned int count = 0;
 unsigned char character;
 if (fm_index->alphabetically_encoded == 1)
  character = c;
 else
  character = get_alphabet_index(fm_index->alphabet,c);
 while(bwt_position%fm_index->sample_SA_size!=0)
 {
  c = fm_index->bwt[bwt_position];
  //printf("bwt position %d znak %c, count %d\n",bwt_position,c,count);
  bwt_position = last_to_first_encoded(c, bwt_position, fm_index->alphabet,fm_index->count_table,fm_index->bwt,fm_index->occurence_table,fm_index->sample_OCC_size);
  count++;
 }
 //printf("bwt position %d znak %c, count %d\n",bwt_position,c,count);
 unsigned int result = fm_index->sampleSA[bwt_position/fm_index->sample_SA_size]+count;
 if (result>fm_index->length)
  result = result - fm_index->length;
 /*if (!(bwt_position)){
  if (count)
    return --count;
  else 
    return count;
}*/
 return result;
}

unsigned int get_SA_value_compressed(unsigned int bwt_position, unsigned int length, unsigned int sample_SA_size, unsigned char *alphabet, unsigned int*sampleSA, 
  unsigned int *count_table,unsigned int block_size, struct compressed_block*block, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman)
{
 //printf("%c znak\n",c);
 unsigned int count = 0;
 unsigned int index = bwt_position/block_size;
 unsigned int position;
 unsigned char c;
 while(bwt_position%sample_SA_size!=0)
 {
  position = bwt_position - index*block_size;
  unsigned char *bwt = decompress_block(block[index].bitvector_length,block[index].bitvector,flag_mtf,flag_runs,flag_huffman,position+1,alphabet,block[index].huffman_tree);
  c = bwt[position];
  free(bwt);
  //printf("bwt position %d znak %c, count %d\n",bwt_position,c,count);
  bwt_position = last_to_first_in_compressed_FMIndex(c,bwt_position,alphabet,count_table,block, block_size, flag_mtf,flag_runs,flag_huffman);
  //printf("**********-nova wt position %d znak %c, count %d\n",bwt_position,c,count);

  index = bwt_position/block_size;
  count++;
 }

 //printf("bwt position %d znak %c, count %d\n",bwt_position,c,count);
 unsigned int result = sampleSA[bwt_position/sample_SA_size]+count;
 if (result>=length)
  result = result - length;
 /*if (!(bwt_position)){
  if (count)
    return --count;
  else 
    return count;
}*/
 return result;
}

unsigned int last_to_first(unsigned char c, unsigned int bwt_position, unsigned char *alphabet,unsigned int *count_table,unsigned char*bwt,unsigned int**occurence_table,unsigned int sample_OCC_size)
{
 unsigned char character = get_alphabet_index(alphabet,c);
 unsigned int last = count_table[character] + count_occ(bwt,occurence_table,bwt_position,c,character,sample_OCC_size);
 //printf("predchodca rotacie %d je %d\n",bwt_position,last);
 return last;
}


unsigned int last_to_first_in_compressed_FMIndex(unsigned char c, unsigned int bwt_position, unsigned char*alphabet, unsigned int* count_table,struct compressed_block *block, unsigned int block_size, 
  unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman)
{
 unsigned int character = get_alphabet_index(alphabet,c);
 unsigned int last = count_table[character] + count_occ_in_compressed_FMIndex(block,block_size,bwt_position,c,flag_mtf, flag_runs, flag_huffman, alphabet);
 //printf("predchodca rotacie %d je %d\n",bwt_position,last);
 return last;
}

unsigned int last_to_first_encoded(unsigned char c, unsigned int bwt_position, unsigned char *alphabet,unsigned int *count_table,unsigned char*bwt,unsigned int**occurence_table,unsigned int sample_OCC_size)
{
 unsigned int last = count_table[c] + count_occ(bwt,occurence_table,bwt_position,c,c,sample_OCC_size);
 //printf("predchodca rotacie %d je %d\n",bwt_position,last);
 return last;
}


unsigned int *create_count_table(unsigned char *s, unsigned int string_length, unsigned char* alphabet)
{
 unsigned int alphabet_size = strlen(alphabet);
 unsigned int *count_table = (unsigned int *) calloc (alphabet_size,sizeof(unsigned int));
 int i = 0;
 int j = 0;
 int k = 0;
 for (i=0;i<string_length;i++)
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

unsigned int **create_occurence_table(unsigned char *s, unsigned int string_length, unsigned char *alphabet, unsigned int sample_size)
{
 unsigned int alphabet_size = strlen(alphabet);
 unsigned int samples_count = (string_length+sample_size)/sample_size;
 unsigned int i;
 unsigned int j;
 unsigned int index = 0;
 unsigned int **occurence_table = (unsigned int **)malloc(alphabet_size*sizeof(unsigned int*));
 for (i=0;i<alphabet_size;i++)
 {
  occurence_table[i] = (unsigned int *)malloc((samples_count+1)*sizeof(unsigned int));
  if (occurence_table[i]==NULL)
  {
   //printf("Error. No memory when allocating occurence table\n");
   exit(1);
  }
 }
 //initialize number of occurences at zero position
 for (i=0;i<alphabet_size;i++)
 {
   occurence_table[i][0] = 0;
 }
 for (i=0;i<string_length;i++)
 {
  index = (i+sample_size-1)/sample_size;
  if (i%sample_size==1)
   for (j=0;j<alphabet_size;j++)
    occurence_table[j][index]=occurence_table[j][index-1];
  for (j=0;j<alphabet_size;j++)
   if (s[i]==alphabet[j])
    occurence_table[j][index]++;     
 } 
 return occurence_table;
}

unsigned int count_occ(unsigned char *s, unsigned int **occurence_table, unsigned int position, unsigned char c, unsigned char character, unsigned int sample_size)
{
 unsigned int count = 0;
 //printf("hladam pre %d, c je %d, char je %d\n",position,c,character);

 unsigned int a = position/sample_size;
 unsigned int bound = position - (a * sample_size);
 //printf("dolne bound je %d, pos je %d\n",bound,position);
 while (bound>0)
 {
  bound--;
  position--;
  //printf("pos %d znak %d %d, akt count je %d\n",position, s[position],c,count);
  if (s[position]==c)
   count++;
 }
 //printf("count je %d, position je %d\n",count,position);
 if (s[position]==c)
   count--;
 //printf("do %d je pocet %c rovny %d, count=%d occ=%d\n",position,c,occurence_table[character][position/sample_size]+count,occurence_table[character][position/sample_size],count);
 //printf("vraciam %d + %d\n",occurence_table[character][position/sample_size],count);
 return occurence_table[character][position/sample_size]+count;
}

unsigned int count_occ_in_block(struct compressed_block *block, unsigned int position,unsigned char c,unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman, unsigned char*alphabet)
{
 unsigned int i;
 unsigned int count = 0;
 unsigned int alphabet_index = get_alphabet_index(alphabet,c);
 unsigned int occurences = block->occurences[alphabet_index];
 //printf("occ %d je %d\n",c,occurences);

 unsigned char *bwt = decompress_block(block->bitvector_length,block->bitvector,flag_mtf,flag_runs,flag_huffman,position,alphabet,block->huffman_tree);
 for (i=0;i<position;i++)
  if (bwt[i] == c)
    count++;
 //printf("\nvraciam %d  + %d\n",occurences,count);
 free(bwt);
 return (occurences+count);
}

unsigned int count_occ_in_compressed_FMIndex(struct compressed_block *block, unsigned int block_size, unsigned int position, unsigned char c, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman, unsigned char*alphabet)
{
 unsigned int index = position/block_size;
 unsigned int remainder = position - index*block_size; //faster than modulo ?
 //printf("pozeram index %d a rem %d,block_size %d pre poziciu %d\n",index,remainder,block_size,position);
 return count_occ_in_block(&block[index],remainder,c,flag_mtf, flag_runs, flag_huffman, alphabet);
}

unsigned int count_occ_in_decompressed_FMIndex(struct compressed_block *block, unsigned char*bwt,unsigned int block_size, unsigned int position, unsigned char c, unsigned char*alphabet)
{
  unsigned int index = position/block_size;
  unsigned int start = index*block_size;
  int count = 0;
  printf("index je %d, start je %d, position %d, c je %c\n",index,start,position,c);
  while (start<position)
  {
    //printf("start=%d, c = %c\n",start,bwt[start]);
    if (bwt[start]==c)
      count++;
    start++;
  }
  if (bwt[position]==c)
      count++;
  unsigned char a = get_alphabet_index(alphabet,c);
  printf("occ %d + count %d\n",block[index].occurences[a],count);
  return (block[index].occurences[a]+count);
}

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
  bwt = (unsigned char *) realloc(bwt,index+cb->block_size);
  strcpy(&bwt[index],decompressed);
  index = index + cb->block_size;
  //free(compressed_FM_index->array_of_blocks[i].occurences);
  free(compressed_FM_index->array_of_blocks[i].bitvector); 

  printf("\n");
  
 }
 //printf("velkost bwt je %d\n",index);
 printf("%s\n",bwt);
 
//printf("Reversing .... \n");
reverseBWT_compressed(bwt,index,compressed_FM_index->end,compressed_FM_index->count_table, alphabet,block_size,compressed_FM_index->array_of_blocks,flag_mtf,flag_runs,flag_huffman);
 // FREE
 return bwt;
}

unsigned char* reverseBWT_compressed(unsigned char*bwt, unsigned int length, unsigned int end, unsigned int* count_table, unsigned char*alphabet, 
  unsigned int block_size, struct compressed_block *block, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman)
{
  unsigned wow = end;
 unsigned int i=0;
 unsigned int j = length-1;
 unsigned char a;
 unsigned char character;
 printf("length je %d\n",length);
 unsigned char *reversed = (unsigned char *)malloc(length);
 if (reversed==NULL)
 {
  printf("Error. No memory when allocating memory when reversing string\n");
  exit(1);
 }
 while(i++!=length)
 {
  a = bwt[end];
  reversed[j] = a;
  j--;
  printf("end=%d,a je %d=%c\n",end,a,a);
  character = get_alphabet_index(alphabet,a);
  printf("pricom %d + %d\n",count_table[character],count_occ_in_decompressed_FMIndex(block,bwt,block_size,end,a,alphabet));
  end = count_table[character]-1 + count_occ_in_decompressed_FMIndex(block,bwt,block_size,end,a,alphabet);
  if (end>length){
    printf("end je %d, length je %d\n",end,length);
    end = end - length;
  }
  //printf("nova pozicia je %d\n",end);
 }
 printf("REVERSED JE\n%s\n",reversed);
}

void alphabet_decode (unsigned char *s, unsigned char * alphabet)
{
 unsigned int string_length = strlen(s);
 unsigned int i;
 for (i=0;i<string_length;i++)
 {
  s[i] = alphabet[s[i]];
 }
}

/*CAGGAGCCGGGTCAAGTCTTGCGCGCGTACCTTCGAAGCGGGCGAACAGGAAGCATGACTTCGGCTGGGCCCCGTCCCCACCGGGCACGCCGTAGAGGGGCAGGGTCATCGGCGCTTGAAAACTGGGGGAA
TGGCCGGCCAGCGTTGCACGACGGGCGCGTTGTGATTAGGCGCCTAGGGGCTGGCCTTGCACCGCGCTCATTCAAGGGGAGACTACAAGACGAGGGCCCCCTCCACACTTGGCTGTAGGAACGTACCAGCTGC
TCTTCCCAAGTCTGGAGGGGGGCCCGACTCGTAACTCCCATTCCGCACGGCCGGAAGCGACCGAGCGTGTCGGCAGCCCCGGCTGACGCATGAAACAACTCGATTGCGGAGACGCCTCAG*/
unsigned char *reverseBWT(unsigned int length, unsigned int end, unsigned char*alphabet,unsigned int*count_table,unsigned char*bwt, unsigned int**occurence_table,unsigned int sample_OCC_size,unsigned char alphabetically_encoded)
{
 int i = 0;
 int j = length -1;
 char a;
 char *reversed = (char *)malloc(sizeof(char)*length);
 if (reversed==NULL)
 {
  printf("Error. No memory when allocating memory when reversing string\n");
  exit(1);
 }
 if (alphabetically_encoded == 0)
 {
  while (i++!=length)
  {
   a = bwt[end];
   reversed[j] = a;
   j--;
   end = last_to_first(a,end,alphabet,count_table,bwt,occurence_table,sample_OCC_size);
  }
 }
 else
 {
  while (i++!=length)
  {
   a = bwt[end];
   reversed[j] = a;
   j--;
   end = last_to_first_encoded(a,end,alphabet,count_table,bwt,occurence_table,sample_OCC_size);
  }
  alphabet_decode(reversed,alphabet);
 }
return reversed;
}


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
 unsigned int sample_OCC_size = fm_index->sample_OCC_size;
 unsigned int length = fm_index->length;
 printf("Printing occurence table:\n");
 for (j=0;j<(length+sample_OCC_size)/sample_OCC_size;j++)
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
 int samples = (fm_index->length+fm_index->sample_SA_size)/fm_index->sample_SA_size;
 int *sample_SA = fm_index->sampleSA;
 printf("Printing sample suffix array:\n");
 for (i = 0; i<samples;i++)
  printf("%d = %d\n",i*fm_index->sample_SA_size,sample_SA[i]);
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
 printf("Length of BWT is: %d\n",fm_index->length);
 printf("Stored BWT is: %s--\n",fm_index->bwt);
 print_bit_vector(fm_index->bwt, fm_index->bitvector_length);
 printf("Reversed   is: %s--\n",reverseBWT(fm_index->length, fm_index->end,fm_index->alphabet,fm_index->count_table,fm_index->bwt, fm_index->occurence_table,fm_index->sample_OCC_size,fm_index->alphabetically_encoded));
 print_count_table(fm_index);
 print_sample_SA(fm_index);
 print_occurence_table(fm_index);
}

unsigned int*search_pattern(struct FMIndex *fm_index, char *pattern)
{
 unsigned int*result=(unsigned int*)malloc(2*sizeof(unsigned int));
 unsigned char*alphabet = fm_index->alphabet;
 unsigned char alphabet_size = strlen(alphabet)-1;
 unsigned int *count_table = fm_index->count_table;
 unsigned int pattern_length = strlen(pattern);
 int j = pattern_length - 1;
 int alphabet_index = get_alphabet_index(alphabet,pattern[j]);
 printf("prvy znak %c index %d\n",pattern[j],alphabet_index);
 result[0] = count_table[alphabet_index];
 if (alphabet_index!=alphabet_size)
  result[1] = count_table[alphabet_index + 1]-1;
 else 
  result[1] = fm_index->length-1;
 printf("rozsah je %d az %d \n",result[0],result[1]);
 while (j>0 && result[0]<result[1])
 {
  j--;
  alphabet_index = get_alphabet_index(alphabet,pattern[j]);
  //printf("%d-ty znak %c index %d\n",j,pattern[j],alphabet_index);
  result[0] = count_table[alphabet_index] + count_occ(fm_index->bwt,fm_index->occurence_table,result[0],alphabet_index,alphabet_index,fm_index->sample_OCC_size);
  result[1] = count_table[alphabet_index] + count_occ(fm_index->bwt,fm_index->occurence_table,result[1],alphabet_index,alphabet_index,fm_index->sample_OCC_size);
 
 /*int test = result[0];
 printf("-----------------------------\n");
 while (test<result[1])
 {
 printf("sa value %d je %d\n",test,get_SA_value(test,fm_index->bwt[result[0]],fm_index));
 test++;
 }
 */
 printf("rozsah je %d az %d \n",result[0],result[1]);

 }
 if (j>0 && result[0]>=result[1])
  result[1]--;

 //printf("vraciam %d %d\n",result[0],result[1]);
 return result;
}

unsigned int*search_pattern_in_compressed_FM_index(struct compressedFMIndex *compressed_fm_index, char *pattern,unsigned char flag_mtf,unsigned char flag_runs, unsigned char flag_huffman)
{
 unsigned int*result=(unsigned int*)malloc(2*sizeof(unsigned int));
 unsigned char*alphabet = compressed_fm_index->alphabet;
 unsigned char alphabet_size = strlen(alphabet)-1;
 unsigned int *count_table = compressed_fm_index->count_table;
 unsigned int pattern_length = strlen(pattern);
 unsigned int block_size = compressed_fm_index->block_size;
 int j = pattern_length - 1;
 int alphabet_index = get_alphabet_index(alphabet,pattern[j]);
 printf("prvy znak %c index %d\n",pattern[j],alphabet_index);
 result[0] = count_table[alphabet_index];
 if (alphabet_index!=alphabet_size)
  result[1] = count_table[alphabet_index + 1]-1;
 else 
 {
  result[1] = compressed_fm_index->genome_length-1;
 }
 printf("\n****rozsah je %d az %d \n",result[0],result[1]);
 while (j>0 && result[0]<result[1])
 {
  j--;
  alphabet_index = get_alphabet_index(alphabet,pattern[j]);
  //printf("%d-ty znak %c index %d\n",j,pattern[j],alphabet_index);
  
  result[0] = count_table[alphabet_index] + count_occ_in_compressed_FMIndex(compressed_fm_index->array_of_blocks,block_size,result[0],pattern[j],flag_mtf,flag_runs,flag_huffman,alphabet);
  result[1] = count_table[alphabet_index] + count_occ_in_compressed_FMIndex(compressed_fm_index->array_of_blocks,block_size,result[1],pattern[j],flag_mtf,flag_runs,flag_huffman,alphabet);
 
 /*int test = result[0];
 printf("-----------------------------\n");
 while (test<result[1])
 {
 printf("sa value %d je %d\n",test,get_SA_value(test,fm_index->bwt[result[0]],fm_index));
 test++;
 }
 */
 printf("rozsah je %d az %d \n",result[0],result[1]);

 }
 if (j>0 && result[0]>=result[1])
  result[1]--;
 
 unsigned int test = result[0];
 while (test<result[1])
 {
  printf("*************************************\n");
 printf("hladam test %d, gl %d\n",test,compressed_fm_index->genome_length);
 printf("sa value %d je %d\n",test,get_SA_value_compressed(test, compressed_fm_index->genome_length,  compressed_fm_index->sample_SA_size, alphabet,compressed_fm_index->sampleSA,count_table,block_size, compressed_fm_index->array_of_blocks, flag_mtf, flag_runs, flag_huffman));
 test++;
 }

 //printf("vraciam %d %d\n",result[0],result[1]);
 return result;
}

unsigned int*approximate_search(int max_error, struct FMIndex *fm_index, char *pattern)
{
 unsigned int pattern_length = strlen(pattern);
 unsigned int div_pattern_length = (pattern_length+1)/(max_error+1);
 unsigned int i,j; 
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
  printf("pozicie su %d az %d\n",result[0],result[1]);
 }
return result;
}

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

unsigned int score(char a, char b)
{
  if (a==b)
    return 1;
  else 
    return -1;
}

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
   printf("%d i je %d, c %d %d\n",nearestMax,i,currentX,currentY);
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
}
