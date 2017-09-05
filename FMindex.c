#include <stdio.h>
#include <stdlib.h>
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
