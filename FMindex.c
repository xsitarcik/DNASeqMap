#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FMindex.h"

int *create_sample_SA(int *suffix_array,int sample_size, int array_size)
{
 int i = 0;
 int sampled_index=0;
 int samples_count = (array_size+sample_size)/sample_size;
 int *sampleSA = (int *) malloc((samples_count)*sizeof(int)); 
 printf("velkost sample pola je: %d",samples_count);
 while (i!=array_size)
 {
  if (i%sample_size == 0)
  {
   sampleSA[sampled_index++]=suffix_array[i];
   printf("uchovam %d hodnotu %d na poziciu %d\n ",i,suffix_array[i],sampled_index);
  }
  i++;
 }
 printf("sapmled je : \n");
 for (i=0;i<samples_count;i++)
  printf("sa[%d]: %d\n",i*sample_size,sampleSA[i]);
 return sampleSA; 
}

//procedure to get SA value from sample suffix array
int get_SA_value(int *sample_SA, int sample_size, int bwt_position)
{
 int count = 0;
 while(bwt_position%sample_size!=0)
 {
  bwt_position = last_to_first(bwt_position);
  count++;
 }
 return (sample_SA[bwt_position%sample_size]+count);
}

int last_to_first(int bwt_position)
{
 return bwt_position;
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
 int **occurence_table = (int **)malloc(alphabet_size*sizeof(int));
 int i;
 int j;
 int index = 0;
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
  printf("i je %d, index je %d\n",i,index);
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
 printf("pos je %d, turn %d\n",position,turn);
 while (position%sample_size!=0)
 {
  printf("position je %d, count je %d \n",position,count);
  if (s[position]==c)
   count = count - turn;
  position = position + turn;
 }
 printf("index je %d\n",position/sample_size);
 return occurence_table[character][position/sample_size]+count;
}
