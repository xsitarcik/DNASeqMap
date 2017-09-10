#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwt.h"
#include "FMindex.h"

int genome_length;
int sample_SA_size=2;
int sample_OCC_size=2;

int main(void)
{
 int i,j;
 int *suffix_array = NULL;
 int *sample_SA = NULL;
 char *s = "AACGAT$";
 char *bwt = NULL;
 struct FMIndex *FM_index;
 char *alphabet = "$ACGT";
 //load genome
 //TODO
 genome_length = strlen(s);
 printf("%d genome", genome_length);

 //create suffix array and bwt of genome
 suffix_array = init_suffix_array(suffix_array,s, genome_length);
 bwt = create_bwt(suffix_array,s,genome_length);
 printf("povodny je %s\n",s);
 printf("bwt je: %s\n",bwt);
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
 FM_index->sample_SA_size = sample_SA_size;
 FM_index->sample_OCC_size = sample_OCC_size;
 printf("do 3. pozicie su A:%d C:%d G:%d T:%d\n",count_occ(bwt,FM_index->occurence_table, 5, 'A', 1, sample_OCC_size),count_occ(bwt,FM_index->occurence_table, 5, 'C', 2, sample_OCC_size), count_occ(bwt,FM_index->occurence_table, 5, 'G', 3, sample_OCC_size), count_occ(bwt,FM_index->occurence_table, 5, 'T', 4, sample_OCC_size));
 printf("pred T je %c\n",bwt[last_to_first('T',4,0,FM_index)]);
 printf("pred A je %c\n",bwt[last_to_first('A',1,6,FM_index)]);
 printf("pred G je %c\n",bwt[last_to_first('G',3,3,FM_index)]);
 printf("pred C je %c\n",bwt[last_to_first('C',2,5,FM_index)]);
 printf("pred A je %c\n",bwt[last_to_first('A',1,4,FM_index)]);
 printf("pred A je %c\n",bwt[last_to_first('A',1,2,FM_index)]);
 printf("pred %% je %c\n",bwt[last_to_first('$',0,1,FM_index)]);
//printf("na pozicii 1 %d\n",get_SA_value(1,'$',0,FM_index));
// printf("na pozicii 2 %d\n",get_SA_value(2,'A',1,FM_index));
// printf("na pozicii 3 %d\n",get_SA_value(3,'G',3,FM_index));
 free(suffix_array); 
 printf("printing sample SA\n");
 for (i=0;i<genome_length/sample_SA_size+1;i++)
  printf("%d: %d\n",i,FM_index->sampleSA[i]);
 print_occurence_table(FM_index->occurence_table,strlen(alphabet),sample_OCC_size,genome_length);
 printf("reversed je :%s\n",reverseBWT(bwt,1,genome_length,alphabet,FM_index));
 free(bwt);
return 0;
}
