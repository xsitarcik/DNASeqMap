#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FMindex.h"

struct FMIndex*build_FM_index(int *suffix_array, int sample_SA_size, int sample_OCC_size, int genome_length, char *bwt, char *alphabet)
{
 struct FMIndex *FM_index = (struct FMIndex*) malloc(sizeof(struct FMIndex));
 FM_index->bwt = bwt;
 FM_index->sample_SA_size = sample_SA_size;
 FM_index->sample_OCC_size = sample_OCC_size;
 FM_index->alphabet = alphabet;
 FM_index->length = genome_length;
 FM_index->end = find_end(suffix_array);
 FM_index->sampleSA = create_sample_SA(suffix_array,sample_SA_size,genome_length);
 free(suffix_array);
 FM_index->count_table = create_count_table(bwt,genome_length,alphabet);
 FM_index->occurence_table = create_occurence_table(bwt,genome_length,alphabet,sample_OCC_size);
 return FM_index;
}

int find_end(int *suffix_array)
{
int i = 0;
while(suffix_array[i]!=0)
 i++;
return i;
}

int *create_sample_SA(int *suffix_array,int sample_size, int array_size)
{
 int i = 0;
 int sampled_index=0;
 int samples_count = (array_size+sample_size)/sample_size;
 int *sampleSA = (int *) malloc((samples_count)*sizeof(int)); 
 if (sampleSA==NULL)
  {
   printf("Error. No memory when allocating sample of suffix array\n");
   exit(1);
  }
 while (i!=array_size)
 {
  if (i%sample_size == 0)
   sampleSA[sampled_index++]=suffix_array[i];
  i++;
 }
 return sampleSA; 
}

//procedure to get SA value from sample suffix array
int get_SA_value(int bwt_position, char c, struct FMIndex *fm_index)
{
 //printf("%c znak\n",c);
 int count = 0;
 int character = get_alphabet_index(fm_index->alphabet,c);
 while(bwt_position%fm_index->sample_SA_size!=0)
 {
  c = fm_index->bwt[bwt_position];
//  printf("bwt position %d znak %c, count %d\n",bwt_position,c,count);
  bwt_position = last_to_first(c, bwt_position, fm_index);
  count++;
 }
// printf("bwt position %d znak %c, count %d\n",bwt_position,c,count);
 if (!(bwt_position)){
  if (count)
    return --count;
  else 
    return count;
}
 return (fm_index->sampleSA[bwt_position/fm_index->sample_SA_size]+count);
}

int last_to_first(char c, int bwt_position, struct FMIndex *fm_index)
{
 int character = get_alphabet_index(fm_index->alphabet,c);
 int last = fm_index->count_table[character] + count_occ(fm_index->bwt,fm_index->occurence_table,bwt_position,c,character,fm_index->sample_OCC_size);
 //printf("predchodca rotacie %d je %d\n",bwt_position,last);
 return last;
}

int *create_count_table(char *s, int string_length, char* alphabet)
{
 int alphabet_size = strlen(alphabet);
 int *count_table = (int *) calloc (alphabet_size,sizeof(int));
 int i = 0;
 int j = 0;
 int k = 0;
 for (i=0;i<string_length;i++)
 {
  for (j=0;j<alphabet_size;j++)
  {
   if (s[i]==alphabet[j])
   {
    count_table[j]++;
    break;
   }
  }
 }
 j = count_table[0];
 count_table[0] = 0;
 for (i=1;i<alphabet_size;i++)
 {
  k = count_table[i];
  count_table[i] = count_table[i-1]+j;
  j = k;
 }
 return count_table;
}

int **create_occurence_table(char *s, int string_length, char *alphabet, int sample_size)
{
 int alphabet_size = strlen(alphabet);
 int samples_count = (string_length+sample_size)/sample_size;
 int i;
 int j;
 int index = 0;
 int **occurence_table = (int **)malloc(alphabet_size*sizeof(int*));
 for (i=0;i<alphabet_size;i++)
 {
  occurence_table[i] = (int *)malloc((samples_count+1)*sizeof(int));
  if (occurence_table[i]==NULL)
  {
   printf("Error. No memory when allocating occurence table\n");
   exit(1);
  }
 }
 //initialize number of occurences at zero position
 for (i=0;i<alphabet_size;i++)
 {
   occurence_table[i][0] = 0;
 }
 for (i=0;i<string_length;i++)
 {
  index = (i+sample_size-1)/sample_size;
  if (i%sample_size==1)
   for (j=0;j<alphabet_size;j++)
    occurence_table[j][index]=occurence_table[j][index-1];
  for (j=0;j<alphabet_size;j++)
   if (s[i]==alphabet[j])
    occurence_table[j][index]++;     
 } 
 return occurence_table;
}

int count_occ(char *s, int **occurence_table, int position, char c, int character, int sample_size)
{
 int count = 0;
 //printf("hladam pre %d\n",position);
 //kde je blizsie
 int a = position/sample_size;
 int bound = position - (a * sample_size);
 //printf("dolne bound je %d, pos je %d\n",bound,position);
 while (bound>0)
 {
  bound--;
  position--;
  //printf("pos %d znak %c %c, akt count je %d\n",position, s[position],c,count);
  if (s[position]==c)
   count++;
 }
 //printf("count je %d, position je %d\n",count,position);
 if (s[position]==c)
   count--;
 //printf("do %d je pocet %c rovny %d, count=%d occ=%d\n",position,c,occurence_table[character][position/sample_size]+count,occurence_table[character][position/sample_size],count);
 //printf("vraciam %d + %d\n",occurence_table[character][position/sample_size],count);
 return occurence_table[character][position/sample_size]+count;
}

char *reverseBWT(struct FMIndex *fm_index)
{
 int i = 0;
 int j = fm_index->length -1;
 int end = fm_index->end;
 char a;
 char *reversed = (char *)malloc(sizeof(char)*fm_index->length);
 if (reversed==NULL)
 {
  printf("Error. No memory when allocating memory when reversing string\n");
  exit(1);
 }
 while (i++!=fm_index->length)
 {
  a = fm_index->bwt[end];
  reversed[j] = a;
  j--;
  end = last_to_first(a,end,fm_index);
 }
return reversed;
}

int get_alphabet_index(char *alphabet, char c)
{
 int i = 0;
 while(c!=alphabet[i])
  i++;
 return i;
}

void print_occurence_table(struct FMIndex *fm_index)
{
 int i,j;
 int alphabet_size = strlen(fm_index->alphabet);
 int sample_OCC_size = fm_index->sample_OCC_size;
 int length = fm_index->length;
 printf("Printing occurence table:\n");
 for (j=0;j<(length+sample_OCC_size)/sample_OCC_size;j++)
 {
  printf("%d: ",j*sample_OCC_size);
  for (i=0;i<alphabet_size;i++)
  {
  printf("%d ",fm_index->occurence_table[i][j]);
  }
  printf("\n");
 }
}

void print_sample_SA(struct FMIndex *fm_index)
{
 int i;
 int samples = (fm_index->length+fm_index->sample_SA_size)/fm_index->sample_SA_size;
 int *sample_SA = fm_index->sampleSA;
 printf("Printing sample suffix array:\n");
 for (i = 0; i<samples;i++)
  printf("%d = %d\n",i*fm_index->sample_SA_size,sample_SA[i]);
}

void print_count_table(struct FMIndex *fm_index)
{
 unsigned int i = 0;
 printf("Printing count table: \n");
 for (i = 0; i<strlen(fm_index->alphabet); i++)
  printf("%4c ",fm_index->alphabet[i]);
 printf("\n");
 for (i = 0; i<strlen(fm_index->alphabet); i++)
  printf("%4d ",fm_index->count_table[i]);
 printf("\n");
}

void print_info_fm_index(struct FMIndex *fm_index)
{
 printf("Length of BWT is: %d\n",fm_index->length);
 printf("Stored BWT is: %s\n",fm_index->bwt);
 printf("Reversed   is: %s\n",reverseBWT(fm_index));
 print_count_table(fm_index);
 print_sample_SA(fm_index);
 print_occurence_table(fm_index);
}

unsigned int*search_pattern(struct FMIndex *fm_index, char *pattern)
{
 unsigned int*result=(unsigned int*)malloc(2*sizeof(unsigned int));
 unsigned char*alphabet = fm_index->alphabet;
 unsigned char alphabet_size = strlen(alphabet)-1;
 int *count_table = fm_index->count_table;
 unsigned int pattern_length = strlen(pattern);
 int j = pattern_length - 1;
 int alphabet_index = get_alphabet_index(alphabet,pattern[j]);
 //printf("prvy znak %c index %d\n",pattern[j],alphabet_index);
 result[0] = count_table[alphabet_index];
 if (alphabet_index!=alphabet_size)
  result[1] = count_table[alphabet_index + 1];
 else 
  result[1] = fm_index->length;
 //printf("rozsah je %d az %d \n",result[0],result[1]);
 
 while (j>0 && result[0]<=result[1])
 {
  j--;
  alphabet_index = get_alphabet_index(alphabet,pattern[j]);
  //printf("%d-ty znak %c index %d\n",j,pattern[j],alphabet_index);
  result[0] = count_table[alphabet_index] + count_occ(fm_index->bwt,fm_index->occurence_table,result[0]-1,pattern[j],alphabet_index,fm_index->sample_OCC_size) + 1;
  result[1] = count_table[alphabet_index] + count_occ(fm_index->bwt,fm_index->occurence_table,result[1],pattern[j],alphabet_index,fm_index->sample_OCC_size);
 //printf("rozsah je %d az %d \n",result[0],result[1]);

 }
 if (result[0])
  result[0]--;
 result[1]--;
 return result;
}

unsigned int*approximate_search(int max_error, struct FMIndex *fm_index, char *pattern)
{
 int pattern_length = strlen(pattern);
 int div_pattern_length = (pattern_length+1)/(max_error+1);
 int i,j; 
 int *result;
 printf("%s dlzka vzoru je %d\n",pattern,pattern_length);
 char *patterns[max_error+1];
 for (i=0;i<max_error;i++)
 {
  patterns[i] = (char *)malloc(sizeof(char)*(div_pattern_length+1)); 
  memcpy( patterns[i], &pattern[i*div_pattern_length], div_pattern_length);
  patterns[i][div_pattern_length]= '\0';
 } 
 j = pattern_length-i*div_pattern_length;
 patterns[i] = (char *)malloc(sizeof(char)*(j+1)); 
 memcpy( patterns[i], &pattern[i*div_pattern_length], j);
 patterns[i][j]= '\0';

 for (i=max_error;i>=0;i--)
 {
  printf("%d. vzor je %s\n",i,patterns[i]);
  result = search_pattern(fm_index,patterns[i]);
  printf("pozicie su %d az %d\n",result[0],result[1]);
 }
return result;
}

int calculate(int i, int j, int k)
{
 int ret_value;
 if (i>j)
  if (i>k)
   ret_value = i;
  else 
   ret_value = k;
 else
  if (j>k) 
   ret_value = j;
  else
   ret_value = k;
 if (ret_value<0)
  return 0;
 else
  return ret_value;
}

int score(char a, char b)
{
  if (a==b)
    return 1;
  else 
    return -1;
}

int get_max_array(unsigned char* array, unsigned int length)
{
 int i, max = 0, maxIndex = 0;
 for (i=0;i<length;i++) 
 {
  if (array[i]>max)
  {
    max = array[i];
    maxIndex = i;
  }
 }
 return maxIndex;
}

unsigned char get_nearest_max(unsigned char i, unsigned char j, unsigned char k)
{
 if (i>=j)
  if (i>=k)
   return 1;
  else 
   return 3;
 else
  if (j>=k) 
   return 2;
  else
   return 3;
}

void align(char *p1, char*p2, int error)
{
 unsigned char **matrix;
 unsigned int corrects;
 unsigned int alignment_length;
 char *alignment1;
 char *alignment2;
 int i,j;
 unsigned int nearestMax;
 unsigned int currentX, currentY;
 unsigned int p1length = strlen(p1), p2length = strlen(p2);
 
 //alignment = (char *) malloc((p1length+error+1)*sizeof(char));
 
 matrix = (unsigned char **) malloc(sizeof(unsigned char*) * (p1length + 1));
 for (i=0;i<=p1length;i++)
  matrix[i] = (unsigned char *) calloc(p2length + 1,sizeof(unsigned char)); 
 for (i=0;i<=error;i++)
  matrix[0][i] = 0;
 for (i=0;i<=error;i++)
  matrix[i][0] = 0;
 
 for (i=1;i<=error;i++){
  for (j=1;j<=i+error;j++){
   matrix[i][j] = calculate(matrix[i-1][j-1] + score(p1[i-1],p2[j-1]),matrix[i-1][j] - 1,matrix[i][j-1] - 1);
   printf(" %d ",matrix[i][j]);
  }
  printf("\n");
}
 for (i=1+error;i<=p1length;i++){
  for (j=i-error;j<=i+error;j++){
    matrix[i][j] = calculate(matrix[i-1][j-1] + score(p1[i-1],p2[j-1]),matrix[i-1][j] - 1,matrix[i][j-1] - 1);
    printf(" %d ",matrix[i][j]);
   }
  printf("\n");
}
 corrects = get_max_array(&matrix[p1length][p1length-error],2*error+1);
 corrects = corrects+p1length-error;
 if (p1length-error>matrix[p1length][corrects])
  printf("No alignment found\n");
 else
 {
   alignment_length = p1length+(p1length - matrix[p1length][corrects])+1;
   alignment1 = (char *) malloc(alignment_length*sizeof(char));
   alignment2 = (char *) malloc(alignment_length*sizeof(char));
   alignment1[alignment_length-1] = '\0';
   alignment2[alignment_length-1] = '\0';
   currentY = corrects;
   currentX = p1length;
   for (i = alignment_length-2;i>0;i--)
   {
    if (currentX==0 || currentY==0)
      break;
    nearestMax = get_nearest_max(matrix[currentX-1][currentY-1],matrix[currentX-1][currentY],matrix[currentX][currentY-1]);
    //if (matrix[currentX][currentY]-1 == matrix[currentX-1][currentY-1]) || (matrix[currentX][currentY]+1 == matrix[currentX-1][currentY-1])
    if (nearestMax==1){
      alignment1[i] = p1[currentX-1];
      alignment2[i] = p2[currentY-1];
      currentX = currentX - 1;
      currentY = currentY - 1; 
    }
    else if (nearestMax==2){
      //printf("nearestMax2 je %d\n",matrix[currentX-1][currentY]);
      if (matrix[currentX-1][currentY]-1 == matrix[currentX][currentY]){
       alignment1[i] = p1[currentX-1];
       alignment2[i] = '-';
       currentX = currentX - 1;
      }
      else
      {
       //printf("cant get there %d %d\n",currentX,currentY);
       alignment1[i] = p1[currentX-1];
       alignment2[i] = p2[currentY-1];
       currentX = currentX - 1;
       currentY = currentY - 1; 
      }
    }
    else{
      //printf("nearestMax3 je %d\n",matrix[currentX][currentY-1]);
      if (matrix[currentX][currentY-1]-1 == matrix[currentX][currentY]){
       alignment1[i] = '-';
       alignment2[i] = p2[currentY-1];
       currentY = currentY - 1;
      }
      else
      {
       //printf("cant get there %d %d\n",currentX,currentY);
       alignment1[i] = p1[currentX-1];
       alignment2[i] = p2[currentY-1];
       currentX = currentX - 1;
       currentY = currentY - 1; 
      }
    }

   }
   printf("%d i je %d, c %d %d\n",nearestMax,i,currentX,currentY);
   if (currentX==1)
   {
    if (currentY==1){
     alignment1[i] = p1[0];
     alignment2[i] = p2[0];
    }
    else{
     alignment1[i] = p1[0];
     alignment2[i] = '-';
    }
   }
   else if (currentY==1)
   {
    alignment1[i] = '-';
    alignment2[i] = p2[0];
   } 
   while (i>=0){
    alignment1[i] = ' ';
    alignment2[i] = ' ';
    i--;
   }
   printf("%s\n",alignment1);
   printf("%s\n",alignment2);
   printf("errors: %d\n",alignment_length - p1length-1);
 }
}
