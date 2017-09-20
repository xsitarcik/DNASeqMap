#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FMindex.h"

struct FMIndex*build_FM_index(int *suffix_array, int sample_SA_size, int sample_OCC_size, int genome_length, char *bwt, char *alphabet)
{
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
 return FM_index;
}

int find_end(int *suffix_array)
{
int i = 0;
while(suffix_array[i]!=0)
 i++;
return i;
}

int *create_sample_SA(int *suffix_array,int sample_size, int array_size)
{
 int i = 0;
 int sampled_index=0;
 int samples_count = (array_size+sample_size)/sample_size;
 int *sampleSA = (int *) malloc((samples_count)*sizeof(int)); 
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
int get_SA_value(int bwt_position, char c, struct FMIndex *fm_index)
{
 int count = 0;
 int character = get_alphabet_index(fm_index->alphabet,c);
 while(bwt_position%fm_index->sample_SA_size!=0)
 {
  bwt_position = last_to_first(c, bwt_position, fm_index);
  count++;
 }
 return (fm_index->sampleSA[bwt_position%fm_index->sample_SA_size]+count);
}

int last_to_first(char c, int bwt_position, struct FMIndex *fm_index)
{
 int character = get_alphabet_index(fm_index->alphabet,c);
 int last = fm_index->count_table[character]-1 + count_occ(fm_index->bwt,fm_index->occurence_table,bwt_position,c,character,fm_index->sample_OCC_size);
 return last;
}

int *create_count_table(char *s, int string_length, char* alphabet)
{
 int alphabet_size = strlen(alphabet);
 int *count_table = (int *) calloc (alphabet_size,sizeof(int));
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

int **create_occurence_table(char *s, int string_length, char *alphabet, int sample_size)
{
 int alphabet_size = strlen(alphabet);
 int samples_count = (string_length+sample_size)/sample_size;
 int i;
 int j;
 int index = 0;
 int **occurence_table = (int **)malloc(alphabet_size*sizeof(int*));
 for (i=0;i<alphabet_size;i++)
 {
  occurence_table[i] = (int *)malloc(samples_count*sizeof(int));
  if (occurence_table[i]==NULL)
  {
   printf("Error. No memory when allocating occurence table\n");
   exit(1);
  }
 }
 //initialize number of occurences at zero position
 for (i=0;i<alphabet_size;i++)
 {
  if (s[0]==alphabet[i])
   occurence_table[i][0] = 1;
  else
   occurence_table[i][0] = 0;
 }
 for (i=1;i<string_length;i++)
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

int count_occ(char *s, int **occurence_table, int position, char c, int character, int sample_size)
{
 int count = 0;
 //kde je blizsie
 int a = position/sample_size;
 int first_bound = position - (a * sample_size);
 int second_bound = ((a+1)*sample_size) - position;
 int turn;
 if (first_bound>second_bound)
  turn = 1;
 else turn = -1;
 while (position%sample_size!=0)
 {
  if (s[position]==c)
   count = count - turn;
  position = position + turn;
 }
 return occurence_table[character][position/sample_size]+count;
}

char *reverseBWT(struct FMIndex *fm_index)
{
 int i = 0;
 int j = fm_index->length -1;
 int end = fm_index->end;
 char a;
 char *reversed = (char *)malloc(sizeof(char)*fm_index->length);
 if (reversed==NULL)
 {
  printf("Error. No memory when allocating memory when reversing string\n");
  exit(1);
 }
 while (i++!=fm_index->length)
 {
  a = fm_index->bwt[end];
  reversed[j] = a;
  j--;
  end = last_to_first(a,end,fm_index);
 }
return reversed;
}

int get_alphabet_index(char *alphabet, char c)
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
 int sample_OCC_size = fm_index->sample_OCC_size;
 int length = fm_index->length;
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
 printf("Stored BWT is: %s\n",fm_index->bwt);
 printf("Reversed   is: %s\n",reverseBWT(fm_index));
 print_count_table(fm_index);
 print_sample_SA(fm_index);
 print_occurence_table(fm_index);
}
