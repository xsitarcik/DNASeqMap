#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwt.h"

int genome_length;

int main(void)
{
 int i,j;
 int *suffix_array = NULL;;
 char *s = "ACGTACACGTA";
 char *bwt = NULL;

 //load genome
 //TODO
 genome_length = strlen(s);
 printf("%d genome", genome_length);

 //create suffix array and bwt of genome
 suffix_array = init_suffix_array(suffix_array,s, genome_length);
 bwt = create_bwt(suffix_array,s,genome_length);

 for(i = 0;i<genome_length;++i)
 {
  printf("SA[%2d] = %2d: ", i, suffix_array[i]);
  for (j = suffix_array[i]; j<genome_length;++j)
  {
   printf("%c ", s[j]);
  }
 printf("\n");
 }
 
 free(suffix_array);
 free(bwt);
return 0;
}
