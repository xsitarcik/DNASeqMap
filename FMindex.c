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

 free(suffix_array);
 return sampleSA; 
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

//aux. procedure for getting index of char in alphabet
unsigned char get_alphabet_index(char *alphabet, unsigned char c)
{
 return (unsigned int)(strchr(alphabet,c) - alphabet);
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

//procedure which breaks input pattern and search each separately in FM index - Wavelet tree
long long int approximate_search_in_FM_index_WT(char *pattern, unsigned int*result)
{
 unsigned int pattern_length = strlen(pattern);
 unsigned int k,temp, length,j;
 unsigned char pattern_half[pattern_length+2];
 unsigned int result_length, current_pattern_start;
 char current_error = 0;
 
 //search last part of the pattern
 char searching_part[pattern_length];
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
       char genome_test[pattern_length+2];
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

//procedure which breaks input pattern and search each separately in FM index - Wavelet tree
long long int approximate_search_in_FM_index_entry(char *pattern, unsigned int*result)
{

 /*clock_t t = clock();*/
 unsigned int pattern_length = strlen(pattern);
 unsigned int k,temp, length,j;
 unsigned char pattern_half[pattern_length+2];
 unsigned int result_length, current_pattern_start;
 char current_error = 0;
 
 //search last part of the pattern
 char searching_part[pattern_length];
 unsigned int div_pattern_length = (pattern_length+1)/(max_error+1);
 unsigned int current_pattern_length = pattern_length;

 int i = max_error;
 current_pattern_start = pattern_length;

 //printf("pocitam %s\n",pattern);fflush(stdout);
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
       while ((k<=max_error) && genome[temp]!=pattern[0])
       {
        temp++;
        k++;
       }
       if (k<=max_error)
       {
        char genome_test[pattern_length+2];
        memcpy(genome_test, &genome[temp], pattern_length+1);
        genome_test[pattern_length+1] = '\0';
        if (align(pattern,genome_test)){
          /*t = clock() - t;
          t_approx_search += t;*/
         return temp;
        }
       }
     }
    }
   }
  }
  i--;
 }
 /* t = clock() - t;
  t_approx_search += t;*/

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

//procedure for aligning two strings up to max error
unsigned char align(char *p1, char*p2)
{
 
 //clock_t t = clock();
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

 //printf("i je %d, pn1length, p2length je %d %d, corrects je %d\n",i,p1length,p2length,corrects);
 for (i=p1length-max_error;i<=p1length;i++)
 {
  if ((i-matrix[i][p2length])<=max_error){
    /*t = clock() - t;
          t_align += t;*/
    //printf("koncim align...\n");fflush(stdout);
    return 1;
  }
 }
  
  for (i=p2length-max_error;i<=p2length;i++)
  {
    if ((i-matrix[p1length][i])<=max_error)
    {
      /*t = clock() - t;
      t_align += t;*/
      //printf("koncim align...\n");fflush(stdout);
      return 1;
    }
  }
  /*t = clock() - t;
  t_align += t;*/
  //printf("koncim align...\n");fflush(stdout);
  return 0;
}

void rebuild_FM_index_into_entries(unsigned int*suffix_array, unsigned char*bwt)
{
 unsigned int i = 0, a;
 unsigned int chars_in_entry = 256;
  
  count_table = create_count_table(bwt);
  unsigned int bitword; unsigned int bit_index;

  unsigned int sample_sa_part = chars_in_entry/sample_SA_size;  //4*4 = 16B
  unsigned int root_in_entry = chars_in_entry/8/sizeof(unsigned int); //128bits = 16B
  unsigned int children_in_entry = root_in_entry + 1; //20B

  unsigned int aux_part = 2; //how many in left and right

  //12Bytes left = 
  //max bitov v unsigned int je 32, log32 +1 = 6bits 6*4 = 24bits = 3B
  //we should need 8Bytes for OCC children 
  unsigned int occ_root_part = 2; //2*4 = 8
  unsigned int total = (genome_length+255) / 256;
  
  entries = (unsigned int*) calloc((total+1)*32,sizeof(unsigned int)); //128/4 = 32
  printf("For entries was totally allocated %u bytes\n", (total+1)*32*4);
  
  //iterators
  unsigned int sa_iter = 0, bit_iter = 0, left_iter = 0,right_iter = 0,root_counter_iter = 0,left_counter_iter = 0, right_counter_iter = 0;
  unsigned int carry = 0;
  unsigned int left_count = 0;
  unsigned int right_count = 0;
  unsigned int temp_count;
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
 //store wt root bitvector, calc how many 1s and 0s
 right_count = 0;
 for (a = 0; a<root_in_entry/2;a++)
 {
  entries[i++] = root->bitvector[bit_iter]>>32; 
  entries[i++] = root->bitvector[bit_iter];
  right_count += __builtin_popcountll(root->bitvector[bit_iter]); //how many 1s
  bit_iter++;
 }
 left_count = chars_in_entry - right_count;
 
 //store children info about how many
 entries[i++] = (left_count+32)/32; //TOTO NIE JE POTREBNE
 
 //store Left wt_child bitvectors
 bitword = 0; bit_index = 0;
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
    entries[i++] = bitword;
    bit_index = 0;
    bitword = 0;
   }
  }
 }

//shift left remaining bits
if (bit_index){
 entries[i++] = bitword << (32-bit_index);
}
else
  entries[i++] = 0;

//store right bitvector child
 bitword = 0; bit_index = 0;
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
    entries[i++] = bitword;
    bit_index = 0;
    bitword = 0;
   }
  }
 }

//shift left remaining bits
if (bit_index) {
 entries[i++] = bitword << (32-bit_index);
 }
 

 //store counters of each ROOT word in a byte
 temp_count = 0;
  for (a = 0; a < 4; a++)
  {
    entries[i] = entries[i] << 8;
    entries[i] += temp_count;
    carry = __builtin_popcount (entries[i+a-18]); 
    temp_count += carry;
  }
  i++;
  for (a = 4; a < 8; a++)
  {
    entries[i] = entries[i] << 8;
    entries[i] += temp_count;
    carry = __builtin_popcount (entries[i+a-19]); 
    temp_count += carry;
  }
  i++;

  //store left counters
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
 
  //store right counters
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

 //storing SA_values
 for (a = 0; a<sample_sa_part;a++)
  entries[i++] = sampleSA[sa_iter++];
 

 }

//handling last part
 right_count = 0;
 for (a = 0; a<root_in_entry/2;a++)
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
 
 //store children info about how many
 entries[i++] = (left_count+32)/32; 
 
 //store wt_children positions
 bitword = 0; bit_index = 0;
 for (a = (i/32)*256; a<genome_length;a++)
 {
  if (bwt[a]=='A' || bwt[a]=='T')
  {
   if (bwt[a]=='T')
   {
    bitword = bitword<<1; bitword = bitword + 1;
    bit_index++;
   }
   else
   {
    bitword = bitword<<1;
    bit_index++;
   }

   if (bit_index==32)
   {
    entries[i++] = bitword;
    bit_index = 0; bitword = 0;
   }
  }
 }
if (bit_index)
 entries[i++] = bitword << (32-bit_index);
else
  entries[i++] = 0;

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
    entries[i++] = bitword;
    bit_index = 0;
    bitword = 0;
   }
  }
 }

if (bit_index) {
 entries[i++] = bitword << (32-bit_index);
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

 //storing SA_values
 for (a = 0; a<sample_sa_part;a++)
  entries[i++] = sampleSA[sa_iter++];

  free(sampleSA);
  
 //building kmers hash table
 total_kmers = 1<<(k_mers_permutation*2); 
 kmers_hash = (unsigned int*) malloc ((total_kmers+1)*sizeof(unsigned int));
 unsigned int *result = (unsigned int *) malloc (2*sizeof(unsigned int));
 char *dna_seq = (char*)malloc(sizeof(char)*(k_mers_permutation+1));
 dna_seq[k_mers_permutation]='\0';

 clock_t t = clock();
 for (i=0;i<total_kmers;i++){
    index_to_dna(i,dna_seq);
    search_pattern_in_FM_index_entry_without_hash(dna_seq,k_mers_permutation,result);
    kmers_hash[i] = result[0];
  }
  kmers_hash[i] = genome_length-1;

  t = clock() - t;
  double time_spent1 = (double)(t) / CLOCKS_PER_SEC;
  printf("It took %lf seconds to build kmer hash map index\n", time_spent1);
  free(result);
  free(dna_seq);
  free(root);
  free(left);
  free(right);
  free(bwt);
}

//procedure for searching pattern in compressed BWT of FM Index
unsigned char search_pattern_in_FM_index_entry_without_hash(char *pattern, unsigned char current_pattern_length, unsigned int*result)
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
   
   unsigned int entry_index1 = result[0]/256; //in each entry up to 256 chars
   unsigned int in_entry_index1 = result[0]%256; //v danom entry
   entry_index1 *=32;

   __builtin_prefetch (&entries[entry_index1], 0, 3);
    __builtin_prefetch (&entries[entry_index1+16], 0, 3);

   unsigned int entry_index2 = result[1]/256; //in each entry up to 256 chars
   unsigned int in_entry_index2 = result[1]%256;
   entry_index2 *=32;

  __builtin_prefetch (&entries[entry_index2], 0, 3);
   __builtin_prefetch (&entries[entry_index2+16], 0, 3);

   result[0] = count_table[alphabet_index] + wt_rank_entry(alphabet_index, entry_index1, in_entry_index1);
   result[1] = count_table[alphabet_index] + wt_rank_entry(alphabet_index, entry_index2, in_entry_index2);
  }

 if (result[0]>=result[1])
 {
  result[1]--;
  return j;
 }
  
 return j;
}

//procedure for searching pattern in compressed BWT of FM Index
unsigned char search_pattern_in_FM_index_entry(char *pattern, unsigned char current_pattern_length, unsigned int*result)
{

  //clock_t t = clock();
 int j = current_pattern_length - k_mers_permutation;
 unsigned char alphabet_index;// = get_alphabet_index(alphabet,pattern[j]);
 unsigned int start,end;
 unsigned int entry_index1,entry_index2,in_entry_index1,in_entry_index2;
 unsigned int index = dna_to_index(&pattern[current_pattern_length - k_mers_permutation]);

 result[0] = kmers_hash[index];
 result[1] = kmers_hash[index+1];

//printf("--------------\n");
 while (result[0]<result[1] && j )
  {
   alphabet_index = get_alphabet_index(alphabet,pattern[--j]);
   
   //printf("res: %d %d\n",result[0], result[1]);
   /*unsigned int entry_index1 = (result[0]>>8) * 32; //in each entry up to 256 chars
   unsigned int in_entry_index1 = result[0]&255; //v danom entry

   unsigned int entry_index2 = (result[1]>>8) * 32; //in each entry up to 256 chars
   unsigned int in_entry_index2 = result[1]&255; //v danom entry*/

   entry_index1 = result[0]/256; //in each entry up to 256 chars
   in_entry_index1 = result[0]%256; //v danom entry
   entry_index1 *=32;

   __builtin_prefetch (&entries[entry_index1], 0, 2);
   __builtin_prefetch (&entries[entry_index1+16], 0, 2);

   entry_index2 = result[1]/256; //in each entry up to 256 chars
   in_entry_index2 = result[1]%256;
   entry_index2 *=32;

  __builtin_prefetch (&entries[entry_index2], 0, 2);
  __builtin_prefetch (&entries[entry_index2+16], 0, 2);

   result[0] = count_table[alphabet_index] + wt_rank_entry(alphabet_index, entry_index1, in_entry_index1);
   result[1] = count_table[alphabet_index] + wt_rank_entry(alphabet_index, entry_index2, in_entry_index2);
  }

 /*t = clock() - t;
 t_access += t;*/
//printf("res: %d %d\n",result[0], result[1]);
 if (result[0]>=result[1])
 {
  //printf("JE TO POTREBNE?\n");
  result[1]--;
  return j;
 }
  
 return j;
}


unsigned int dna_to_index(char *s){
 
 unsigned int i = 0;
 unsigned int result = 0;
 //printf("string je %s\n",s);
//4^10 - 4MB
 for (i=0;i<k_mers_permutation;i++) {
  result = result << 2;
  result = result + get_alphabet_index(alphabet,s[i]);
  /*switch(s[i])
  {
    case 'C': 
      result = result + 1; //01
      break;
    case 'G': 
      result = result + 2; //10
      break;
    case 'T': 
      result = result + 3; //11
      break;
  }*/
 }
 return result;
}

void index_to_dna(unsigned int index, char *dna){
 
 unsigned int i,c;
//4^10 - 4MB
 for (i=0;i<k_mers_permutation;i++) {
  c = index & 0x03;
  index = index >>2;
  switch(c)
  {
    case 0 :
      dna[k_mers_permutation-i-1] = 'A'; //00
      break;
    case 1:
      dna[k_mers_permutation-i-1] = 'C'; //01
      break;
    case 2:
      dna[k_mers_permutation-i-1] = 'G'; //10
      break;
    case 3:
      dna[k_mers_permutation-i-1] = 'T'; //11
      break;
  }
 }
}

unsigned int wt_rank_entry(unsigned char c, unsigned int entry_index, unsigned int in_entry_index)
{

  //clock_t t = clock();

 unsigned int i;
 unsigned int count = 0;
 unsigned int result = 0;
 unsigned int a;
 unsigned int remainder;
 

 i = in_entry_index/32;
 remainder = in_entry_index%32;
 if (i>3){
    count = entries[entry_index+19]>>((7-i)*8) & 0xff;
  }
 else
  {
    count = entries[entry_index+18]>>((3-i)*8) & 0xff; 
  }

 if (remainder)
   count = count + __builtin_popcount(entries[entry_index+i]>>(32-remainder)); 

 
 if (c == 0 || c == 3)
 {
  //bolo nutne spocitat nuly, preto invertujem
  count = in_entry_index - count; //VIEM KOLKO BOLO NUL DO POZICIE POS V ROOT

  //spocitam 1 do pozicie count
  remainder = count;
  for (i=0;i<count/32;i++)
    {
     result += __builtin_popcount(entries[entry_index+9+i]);
     remainder -= 32;
    }
  
   if (remainder){
      result += __builtin_popcount(entries[entry_index+9+i]>>(32-remainder));
  }
  //VIEM POCET 1 OD ZACIATKU ENTRY PO POS
  // v entry indexe je pocet 1
  //vysledok je = count
  /* t = clock() - t;
  t_rank += t;*/

  if (c)
   return (result + entries[entry_index + 21]);
  else
   return (count - result + entries[entry_index + 20]);
 }

 else 
 {
  a = entries[entry_index + 8] + entry_index + 9; //kolko bolo v left?
  remainder = count;
  for (i=0;i<count/32;i++)
  {
   result = result + __builtin_popcount(entries[a+i]);
   remainder -= 32;
  }

  if (remainder){
   result = result + __builtin_popcount(entries[a+i]>>(32-remainder));
 }
 /*t = clock() - t;
 t_rank += t;*/


 if (c==2)
   return (count - result + entries[entry_index + 22]);
  else
   return (result + entries[entry_index + 23]);
 }
}

unsigned int wt_access_rank_entry(unsigned int entry_index, unsigned int in_entry_index)
{

 unsigned int i;
 unsigned int count = 0;
 unsigned int result = 0;
 unsigned int a = in_entry_index;
 unsigned int remainder = in_entry_index;
 
 i = in_entry_index/32;
 remainder = in_entry_index%32;
 if (i>3){
    count = entries[entry_index+19]>>((7-i)*8) & 0xff;
  }
 else{
    count = entries[entry_index+18]>>((3-i)*8) & 0xff; 
  }

 if (remainder)
   count = count + __builtin_popcount(entries[entry_index+i]>>(32-remainder)); 

  if (get_bit_4byte(entries[entry_index+i],remainder))
  {

    //printf("pocitam 1, entry_index %d\n",entry_index);
    //aky znak na count pozicii?
    a = entries[entry_index + 8] + entry_index + 9; //kolko bolo v left?
    remainder = count;
    for (i=0;i<count/32;i++)
    {
     result = result + __builtin_popcount(entries[a+i]);
    remainder -= 32;
    }

    if (remainder)
      result = result + __builtin_popcount(entries[a+i]>>(32-remainder));
    
    //remainder = count%32;

    //printf("pre count %u, pristupujem na poziciu %u %u\n",count,count/sample_SA_size,remainder);
    //printf("entries[%d] je %u\n",entry_index + 9 + a + count/32, entries[entry_index + 9 + a + count/32]);
    if (get_bit_4byte(entries[a + count/32], remainder))
      return (result + entries[entry_index + 23])+count_table[1]; //C 11
    else
      return (count - result + entries[entry_index + 22])+count_table[2]; //G 10
  }

 else
 {
  //bolo nutne spocitat nuly, preto invertujem
  count = in_entry_index - count;

  //spocitam 1 do pozicie count
  remainder = count;
  for (i=0;i<count/32;i++)
    {
     result += __builtin_popcount(entries[entry_index+9+i]);
     remainder -= 32;
    }
  
   if (remainder)
    result += __builtin_popcount(entries[entry_index+9+i]>>(32-remainder));

  if (get_bit_4byte(entries[entry_index + 9 + count / 32], remainder))
   return (result + entries[entry_index + 21])+count_table[3];//T - 01
  else
   return (count - result + entries[entry_index + 20])+count_table[0]; //A 00
 }
}

unsigned int get_SA_value_entry(unsigned int bwt_position)
{
 unsigned int count = 0;
 unsigned int entry_index; //in each entry up to 256 chars
 unsigned short int in_entry_index; //v danom entry
 //unsigned int character;
 //unsigned char c;

 //clock_t t = clock();

   //TOTO ZEFEKTIVNIT
 while(bwt_position%sample_SA_size!=0)
 {
  entry_index = (bwt_position/256) * 32;
  in_entry_index = bwt_position%256;

  __builtin_prefetch (&entries[entry_index], 0, 0);
  __builtin_prefetch (&entries[entry_index+16], 0, 0);
  bwt_position = wt_access_rank_entry(entry_index, in_entry_index);
  /*character = wt_access_entry(entry_index, in_entry_index);
  bwt_position = count_table[character] + wt_rank_entry(character, entry_index, in_entry_index);*/
  count++;
 }
 
 entry_index = (bwt_position/256) * 32;
 in_entry_index = bwt_position%256;
 count += entries[entry_index + 24 + in_entry_index/32];
 
 if (count>=genome_length)
  count = count - genome_length;


/*t = clock() - t;
t_get_sa_value += t;*/

return count;
}



//function which returns character on input position
unsigned char wt_access_entry(unsigned int entry_index, unsigned short int in_entry_index)
{
  
 unsigned char i;
 unsigned char count = 0;
 unsigned int a = entry_index;
 unsigned char remainder = in_entry_index;
 //printf("indexy su %d %d \n",entry_index,in_entry_index);
  //spocitaju sa 1ky od zaciatku entry po dany in_entry_index

  for (i=0;i<in_entry_index/32;i++,a++)
  {
    count = count + __builtin_popcount(entries[a]);
    //printf("count je %d\n",count);
    remainder = remainder - 32;
  }
  if (remainder){
    //printf("hm remmainder %d\n",remainder);
   count = count + __builtin_popcount(entries[a]>>(32-remainder));
  }
  
  //printf("wut count je %d, i je%d\n",count,i);
  if (get_bit_4byte(entries[a],remainder))
  {

    //printf("pocitam 1, entry_index %d\n",entry_index);
    //aky znak na count pozicii?
    a = entries[entry_index+8];
    
    //remainder = count%32;

    //printf("pre count %u, pristupujem na poziciu %u %u\n",count,count/sample_SA_size,remainder);
    //printf("entries[%d] je %u\n",entry_index + 9 + a + count/32, entries[entry_index + 9 + a + count/32]);
    if (get_bit_4byte(entries[entry_index + 9 + a + count/32], count%32))
      return 1;
    else
      return 2;
  }

  else
  {
    
    //pre nulu je potrebne spocitat nuly, preto invertujem
    count = in_entry_index - count;

    if (get_bit_4byte(entries[entry_index + 9 + count / 32], count % 32))
      return 3;
    else
      return 0;
  }
}

