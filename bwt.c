#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <divsufsort.h>

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

int max(int number1, int number2)
{
 if (number1 < number2)
  return number2;
 else
  return number1;
}

//helper function used when sorting rotations
int compare_rotations(char *s, int start1, int start2)
{
 int rotation_break = max(start1,start2);
 int i = rotation_break;
 int origin1 = start1;
 int origin2 = start2;
 printf("i je: %d, znaky su %d %c a %d %c\n",i,start1,s[start1],start2,s[start2]);
 while (i!=genome_length && s[start1]==s[start2])
 {
  i++;
  start1++;
  start2++;
 }
  printf("porovnavam %d %c, %d %c",start1,s[start1],start2,s[start2]);
 printf("i je: %d, znaky su %d %c a %d %c\n",i,start1,s[start1],start2,s[start2]);
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
 printf("i je: %d, znaky su %d %c a %d %c\n",i,start1,s[start1],start2,s[start2]);
 if (s[start1]<s[start2])
  return origin1;
 else return origin2;
}

char *create_bwt(int *suffix_array, char *s)
{
 int i;
 char *bwt = (char *) malloc (genome_length * sizeof(char));
 printf("porovnavam pozicie 1 a 3: %d\n",compare_rotations(s,1,3));
 for (i = 0; i<genome_length; i++)
  printf("%d %c\n",suffix_array[i],s[i]);
 return bwt;
}

int main(void)
{
 int i,j;
 int *suffix_array = NULL;;
 char *s = "ACGTAGGCGTCCA";
 char *bwt = NULL;
 
 //load genome
 genome_length = strlen(s);

 //create suffix array of genome
 suffix_array = create_suffix_array(suffix_array,s); 
 bwt = create_bwt(suffix_array,s);

 for(i = 0;i<genome_length;++i)
 {
  printf("SA[%2d] = %2d: ", i, suffix_array[i]);
  for (j = 0; j<genome_length;++j)
  {
   printf("%c", s[j]);
  }
 }

 free(suffix_array);
return 0;
}
