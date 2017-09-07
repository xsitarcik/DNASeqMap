#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwt.h"
#include "FMindex.h"

int genome_length;
int sample_SA_size=3;
int sample_OCC_size=4;

int main(void)
{
 int i,j;
 int *suffix_array = NULL;
 int *sample_SA = NULL;
 char *s = "ACGTACACGTA";
 char *bwt = NULL;
 struct FMIndex *FM_index;
 char *alphabet = "ACGT";
 //load genome
 //TODO
 genome_length = strlen(s);
 printf("%d genome", genome_length);

 //create suffix array and bwt of genome
 suffix_array = init_suffix_array(suffix_array,s, genome_length);
 bwt = create_bwt(suffix_array,s,genome_length);
 printf("wtf");
 for(i = 0;i<genome_length;++i)
 {
  printf("SA[%2d] = %2d: ", i, suffix_array[i]);
  for (j = suffix_array[i]; j<genome_length;++j)
  {
   printf("%c ", s[j]);
  }
 printf("\n");
 }
 FM_index = (struct FMIndex*) malloc(sizeof(struct FMIndex));
 FM_index->bwt = bwt;
 FM_index->sampleSA = create_sample_SA(suffix_array,sample_SA_size,genome_length);
 FM_index->count_table = create_count_table(bwt,genome_length,alphabet);
 FM_index->occurence_table = create_occurence_table(bwt,genome_length,alphabet,sample_OCC_size);
 printf("do 3. pozicie su A:%d C:%d G:%d T:%d\n",count_occ(bwt,FM_index->occurence_table, 10, 'A', 0, sample_OCC_size),count_occ(bwt,FM_index->occurence_table, 10, 'C', 1, sample_OCC_size), count_occ(bwt,FM_index->occurence_table, 10, 'G', 2, sample_OCC_size), count_occ(bwt,FM_index->occurence_table, 10, 'T', 3, sample_OCC_size));
 //printf("na pozicii %d\n",get_SA_value(FM_index->sampleSA,sample_SA_size,1));
 //printf("na pozicii %d\n",get_SA_value(FM_index->sampleSA,sample_SA_size,2));
 //printf("na pozicii %d\n",get_SA_value(FM_index->sampleSA,sample_SA_size,3));
 free(suffix_array); 
 for (i=0;i<genome_length/sample_SA_size+1;i++)
  printf("%d: %d\n",i,FM_index->sampleSA[i]);
 free(bwt);
return 0;
}
