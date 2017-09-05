#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwt.h"
#include "FMindex.h"

int genome_length;
int sample_SA_size=3;

int main(void)
{
 int i,j;
 int *suffix_array = NULL;
 int *sample_SA = NULL;
 char *s = "ACGTACACGTA";
 char *bwt = NULL;

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
 //sampe
 sample_SA = create_sample_SA(suffix_array,sample_SA_size,genome_length);
 printf("sapmled je : \n");
 for (i=0;i<genome_length/sample_SA_size+1;i++)
  printf("sa[%d]: %d\n",i*sample_SA_size,sample_SA[i]);
 free(suffix_array);
 free(bwt);
return 0;
}
