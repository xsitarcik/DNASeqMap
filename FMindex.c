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
 printf("samples_count SA je %d\n",samples_count);
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
int get_SA_value(int bwt_position, char c, int character, struct FMIndex *fm_index)
{
 int count = 0;
 while(bwt_position%fm_index->sample_SA_size!=0)
 {
  bwt_position = last_to_first(c, character, bwt_position, fm_index);
  count++;
 }
 return (fm_index->sampleSA[bwt_position%fm_index->sample_SA_size]+count);
}

int last_to_first(char c, int character, int bwt_position, struct FMIndex *fm_index)
{
 printf("mofo\n");
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
 int **occurence_table = (int **)malloc(alphabet_size*sizeof(int*));
 int i;
 int j;
 int index = 0;
printf("s je %s\n",s);
 printf("alphabet size je %d, string length je %d \n",alphabet_size, string_length);
 printf("samples_count je %d\n", samples_count);
 printf("sample_size je %d\n",sample_size);
 for (i=0;i<alphabet_size;i++)
 {
  occurence_table[i] = (int *)malloc(samples_count*sizeof(int));
  if (occurence_table[i]==NULL)
  {
   printf("Error. No memory when allocating occurence table\n");
   exit(1);
  }
 }
 printf("idk\n");
 //initialize number of occurences at zero position
 for (i=0;i<alphabet_size;i++)
 {
  printf("i je %d\n",i);
  if (s[0]==alphabet[i])
   occurence_table[i][0] = 1;
  else
   occurence_table[i][0] = 0;
  printf("done\n");
 }
 printf("wtf ?\n");
 for (i=1;i<string_length;i++)
 {
  index = (i+sample_size-1)/sample_size;
  printf("i je %d a index je %d\n",i,index);
  if (i%sample_size==1)
   for (j=0;j<alphabet_size;j++)
    occurence_table[j][index]=occurence_table[j][index-1];
  printf("%d %d %d %d %d\n",occurence_table[0][index],occurence_table[1][index],occurence_table[2][index],occurence_table[3][index],occurence_table[4][index]);
  for (j=0;j<alphabet_size;j++)
   if (s[i]==alphabet[j])
    occurence_table[j][index]++;     
 } 
 return occurence_table;
}

int count_occ(char *s, int **occurence_table, int position, char c, int character, int sample_size)
{
 printf("position %d, c %c, character %d\n",position,c,character);
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
 printf("vraciam pre poziciu %d a znak %c hodnoty %d\n",position,c,character);
 printf("hodnotu beriem z %d a count je %d\n",position/sample_size,count);
 printf("vysledok je %d\n",occurence_table[character][position/sample_size]+count);
 return occurence_table[character][position/sample_size]+count;
}
char *reverseBWT(char *bwt,int end,int length, char *alphabet, struct FMIndex *fm_index)
{
 int i = 0;
 int j = length -1;
 char *reversed = (char *)malloc(sizeof(char)*length);
 while (i++!=length)
 {
  printf("i je %d , je %d\n",i,j);
  reversed[j]=bwt[end];
  printf("na poziciu %d ukladam %c\n",j,bwt[end]);
  j--;
  end = last_to_first(bwt[end],get_alphabet_index(alphabet,bwt[end]),end,fm_index);
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
void print_occurence_table(int **occurence_table, int alphabet_size, int sample_OCC_size, int length)
{
 int i,j;
 printf("Printing occurence table:\n");
 for (j=0;j<(length+sample_OCC_size)/sample_OCC_size;j++)
 {
  for (i=0;i<alphabet_size;i++)
  {
  printf("%d ",occurence_table[i][j]);
  }
  printf("\n");
 }
}

