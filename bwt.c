#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int genome_length;

int *create_suffix_array(int *suffix_array,char *s)
{
 int i;
 //initialize suffix array
 suffix_array = (int *) malloc(genome_length*sizeof(int));
 if (suffix_array == NULL)
  printf("Error when allocating memory for suffix array\n");
 else
 {
  for (i=0;i<genome_length;i++)
   suffix_array[i] = i;
 }
 
 printf("s je: %s\n",s);
 for (i=0;i<genome_length;i++)
  {
   printf("%d - %c\n",suffix_array[i],s[i]);
  }
 return suffix_array;
}

char *create_bwt(int *suffix_array, char *s)
{
 int i;
 char *bwt = (char *) malloc (genome_length * sizeof(char));
 for (i = 0; i<genome_length; i++)
  printf("%d %c\n",suffix_array[i],s[i]);
 return bwt;
}

int main(void)
{
 int i;
 int *suffix_array = NULL;;
 char *s = "ACGTAGGCGTCCA$";
 char *bwt = NULL;
 
 //load genome
 genome_length = strlen(s);

 //create suffix array of genome
 suffix_array = create_suffix_array(suffix_array,s); 
 bwt = create_bwt(suffix_array,s);
return 0;
}
