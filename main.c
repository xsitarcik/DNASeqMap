#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwt.h"
#include "FMindex.h"
#include "compression.h"

int genome_length;
int sample_OCC_size = 4;
int sample_SA_size = 4;

int main(void)
{
 int i,j;
 int *suffix_array = NULL;
 int *sample_SA = NULL;
 char *s = "AACGATAACGATAACGATTGACAGTA$";
 char *bwt = NULL;
 struct FMIndex *FM_index;
 char *alphabet = "$ACGT";
 struct symbol_table *front = NULL;
 //load genome
 //appendovanie $
 //abeceda
 //TODO
 genome_length = strlen(s);
 //count if memory is sufficient
 // for construction of SA we need cca 5n bytes 
 // then for bwt we need another n bytes
 // so together 6n bytes at least, 
 // if its not sufficient, break origin string TO DO

 //create suffix array and bwt of genome
 suffix_array = init_suffix_array(suffix_array,s, genome_length);
 bwt = create_bwt(suffix_array,s,genome_length);
 //free(genome);
 
 //count if sample sizes are high enough so occ table can fit in memory TO DO
 //build FMIndex and free suffix_array AND GENOME (to do)
 FM_index = build_FM_index(suffix_array,sample_SA_size, sample_OCC_size, genome_length,bwt,alphabet);
 //free(suffix_array); 
 
 //printing
 printf("Org string is: %s\n",s);
 print_info_fm_index(FM_index);
 
 front = build_symbol_table(alphabet);
 move_to_front_encode(front,FM_index->bwt);
//app matching:
// read sequence R one by one
// according to error E
// parse sequence R to E+1 sequences r1... non-overlapping
// each r1... sequence locate with help of FM-INDEX

 free(bwt);

return 0;
}
