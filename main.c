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

 //create suffix array and bwt of genome
 suffix_array = init_suffix_array(suffix_array,s, genome_length);
 bwt = create_bwt(suffix_array,s,genome_length);
 
 //build FMIndex and free suffix_array AND GENOME (to do)
 FM_index = build_FM_index(suffix_array,sample_SA_size, sample_OCC_size, genome_length,bwt,alphabet);
 free(suffix_array); 
 
 //printing
 printf("Original je: %s\n",s);
 printf("BWT je     : %s\n",FM_index->bwt);
 printf("Reversed je: %s\n",reverseBWT(1,FM_index));
 print_occurence_table(FM_index->occurence_table,strlen(alphabet),sample_OCC_size,genome_length);
 print_sample_SA(FM_index);
 free(bwt);
return 0;
}
