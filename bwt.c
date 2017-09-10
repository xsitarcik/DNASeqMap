#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bwt.h"

int *init_suffix_array(int *suffix_array,char *s,int genome_length)
{
 int i;
 //initialize suffix array
 suffix_array = (int *) malloc(genome_length*sizeof(int));
 if (suffix_array == NULL)
 {
  printf("Error when allocating memory for suffix array\n");
  exit(1);
 }
 else
 {
  for (i=0;i<genome_length;i++)
   suffix_array[i] = i;
 }
 return suffix_array;
}

int max(int number1, int number2)
{
 if (number1 < number2)
  return number2;
 else
  return number1;
}

//helper function used when sorting rotations
int compare_rotations(char *s, int start1, int start2, int genome_length)
{
 int rotation_break = max(start1,start2);
 int i = rotation_break;
 int origin1 = start1;
 int origin2 = start2;
 while (i!=genome_length && s[start1]==s[start2])
 {
  i++;
  start1++;
  start2++;
 }
 if (s[start1]==s[start2])
 {
  i = 0;
  while (i!=rotation_break && s[start1]==s[start2])
  {
   i++;
   start1++;
   start2++;
  }
 }
 if (s[start1]<s[start2])
  return origin1;
 else 
  return origin2;
}

int *insertion_sort_array(int *suffix_array, char *s, int genome_length)
{
 int i = 0;
 int j;
 int swap_temp;
 while (i != genome_length)
 {
  j = i;
  while (j > 0 && compare_rotations(s,suffix_array[j-1],suffix_array[j],genome_length)==suffix_array[j])
  {
   swap_temp = suffix_array[j-1];
   suffix_array[j-1] = suffix_array[j];
   suffix_array[j] = swap_temp;
   j--;
  }
  i++;
 }
 return suffix_array;
}

char *create_bwt(int *suffix_array, char *s, int genome_length)
{
 int i;
 char *bwt = (char *) malloc (genome_length * sizeof(char));
 suffix_array = insertion_sort_array(suffix_array,s,genome_length); 
 for(i=0;i<genome_length;i++)
 {
  if (suffix_array[i]==0)
   bwt[i] = s[genome_length-1];
  else bwt[i] = s[suffix_array[i]-1];
 }
 return bwt;
}
