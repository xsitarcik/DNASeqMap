#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bwt.h"

//auxiliary method for loading genome from file
unsigned char*load_genome_from_file(unsigned char*file,unsigned int*genome_length)
{
unsigned char * buffer;
unsigned int length;
FILE * f = fopen (file, "r");

if (f)
{
  fseek (f, 0, SEEK_END);
  length = ftell (f);
  fseek (f, 0, SEEK_SET);
  buffer = malloc (length);
  if (buffer)
  {
    fread (buffer, 1, length, f);
  }
  fclose (f);
}
 buffer[length]='\0';
 *genome_length = length;
 return buffer;
}

//procedura na zotriedenie a zlucenie casti
void topDownMerge(char *s,unsigned int *positions, unsigned int begin, unsigned int mid, unsigned int end, unsigned int *pomocnePole, unsigned int length){
  unsigned int i = begin,j = mid,k;
  for (k = begin; k < end; k++) {

    if(i < mid && (j >= end || (porovnaj(s,positions[i],positions[j],length)))) {
      pomocnePole[k] = positions[i];
            i = i + 1;
        } else {
            pomocnePole[k] = positions[j];
            j = j + 1;    
        }
    }
}

//procedura na rozdelenie na dve casti, ktore sa zotriedia a zlucia
//vyuzitie openMP, vytvorenie paralelnej sekcie pre kazde rozdelenie
void topDownSplitMerge (char *s,unsigned int *positions, unsigned int begin, unsigned int end, unsigned int *pomocnePole, unsigned int length){
  unsigned int mid;
  if ((end - begin) < 2)
    return;
  else {
  mid = (end + begin) / 2;
  #pragma omp parallel sections
        {
        #pragma omp section
  topDownSplitMerge(s,positions, begin, mid, pomocnePole,length);
    #pragma omp section
  topDownSplitMerge(s,positions, mid, end, pomocnePole,length);  
  }
  topDownMerge(s,positions, begin, mid, end, pomocnePole,length);  
  copyArr(pomocnePole,begin,end,positions);
  }
}

//procedura na prekopirovania hodnot pola do druheho pola
void copyArr(unsigned int src[], unsigned int begin, unsigned int end, unsigned int dest[]){
    long long int i;
  #pragma omp parallel for
  for(i = begin; i < end; i++)
        dest[i] = src[i];
}

//procedura na triedenie zlucovanim
void mergeSort(char*s,unsigned int *positions,unsigned int size){
  unsigned int *pomocnePole = (unsigned int*)malloc(size*sizeof(unsigned int));
  omp_set_num_threads(2);
  topDownSplitMerge(s,positions, 0, size, pomocnePole,size);
  free(pomocnePole);
}

unsigned int porovnaj(char *s,unsigned int pos1, unsigned int pos2,  unsigned int length){
  unsigned int i,a,b,e=0;
  a = pos1;
  b = pos2;
  //zistenie ktora pozicia rotacie je blizsia ku koncu retazca
  if ((length-pos1)>(length-pos2)){
    //porovnavaju sa znaky po koniec retazca
    for (i=0;i<length-pos2;i++){
      if (s[a]!=s[b]){
        e = 1;
        break;
      }
      //ak su znaky na danych poziciach rovnake, inkrementujeme pozicie
      a++;
      b++;
    }
    //ak porovnavana pozicia dosiahne koniec retazca
    if (e==0){
      b = 0; //nastavi sa na nulu
      for (i=0;i<pos1;i++){
        if (s[a]!=s[b])
          break;
        a++;
        //ak aj druha porovnavana pozicia dosiahne koniec retazca, nastavi sa na nulu
        if (a==length) 
          a = 0;
        b++;
      }
    }
  }
  //rovnaky postup
  else {
    for (i=0;i<length-pos1;i++){
      if (s[a]!=s[b]){
        e = 1;
        break;
      }
      a++;
      b++;
    }
    if (e==0){
      a = 0;
      for (i=0;i<pos2;i++){
        if (s[a]!=s[b])
          break;
        a++;
        b++;
        if (b==length)
          b = 0;
      }
    }
  }
  
  //vrati 0 ak rotacia zacinajuca na pos1 je lexikograficky vacsia
  //vrati 1 ak rotacia zacinajuca na pos1 je lexikograficky mensia alebo rovna
  if (s[a]>s[b])
    return 0;
  else {
    return 1;
  }
}

//auxiliary method for initiazing array of indexes, which would represent starting positions of rotations
unsigned int *init_suffix_array(unsigned int *suffix_array,unsigned char *s,unsigned int genome_length)
{
 unsigned int i;
 //initialize suffix array
 suffix_array = (unsigned int *) malloc(genome_length*sizeof(unsigned int));
 if (suffix_array == NULL)
 {
  printf("Error when allocating memory for suffix array\n");
  exit(1);
 }
 else
 {
  for (i=0;i<genome_length;i++)
   suffix_array[i] = i;
 }
 return suffix_array;
}

//auxiliary method for getting maximum of 2 numbers
unsigned int max(unsigned int number1,unsigned int number2)
{
 if (number1 < number2)
  return number2;
 else
  return number1;
}

//auxiliary method for getting minimum of 2 numbers
unsigned int min(unsigned int number1,unsigned int number2)
{
 if (number1 > number2)
  return number2;
 else
  return number1;
}

//method for inplace reversing string
void reverse_string(unsigned char *str)
{
    //handle null case
    if (str == 0)
    {
        return;
    }

    //handle empty case
    if (*str == 0)
    {
        return;
    }

    //get start and end pointers
    unsigned char *start = str;
    unsigned char *end = start + strlen(str) - 1;
    unsigned char temp;

    //reverse chars in cycle from start and end
    while (end > start)
    {
        //swap chars
        temp = *start;
        *start = *end;
        *end = temp;

        //move in indexes
        ++start;
        --end;
    }
}

//auxiliary function used when sorting rotations of 2 input starting positions
unsigned int compare_rotations(unsigned char *s,unsigned int start1,unsigned int start2,unsigned int genome_length)
{
 unsigned int rotation_break1 = max(start1,start2) - min(start1,start2);
 unsigned int rotation_break2 = min(start1,start2);
 unsigned int real_length = genome_length - 1;
 unsigned int origin1 = start1;
 unsigned int origin2 = start2;
 unsigned int i = max(start1,start2);

 //compare rotations char by char until end of string
 while (i!=real_length && s[start1]==s[start2])
 {
  i++;
  start1++;
  start2++;
 }

 //if index of one rotation reached end of string, reset index of that rotation to 0
 if (s[start1]==s[start2])
 {
  i = 0;
  if (start1>start2)
  {
   start1 = 0;
   start2++;
  }
  else 
  {
   start1++;
   start2 = 0;
  }

  //compare rotations char by char until second rotation reaches end of string
  while (i!=rotation_break1 && s[start1]==s[start2])
  {
   i++;
   start1++;
   start2++;
  }

  //reset index of second rotation to 0
  if (i == rotation_break1)
  {
   if (start1>start2)
    start1 = 0;
   else 
    start2 = 0;
   i = 0;

   //compare rotations char by char until reached start of earlier rotation
   while (i!=rotation_break2 && s[start1]==s[start2])
   {
    i++;
    start1++;
    start2++;
   }
   if (s[start1]<s[start2])
    return origin1;
   else 
    return origin2;
  }
  else {
   if (s[start1]<s[start2])
    return origin1;
   else 
    return origin2;
  }
 }
 else 
  {
   if (s[start1]<s[start2])
    return origin1;
   else 
    return origin2;
  }
}

//sorting array of rotations with insertion sort
unsigned int *insertion_sort_array(unsigned int *suffix_array, unsigned char *s, unsigned int genome_length)
{
 unsigned int i = 0;
 unsigned int j;
 unsigned int swap_temp;
 while (i != genome_length)
 {
  j = i;
  while (j > 0 && compare_rotations(s,suffix_array[j-1],suffix_array[j],genome_length)==suffix_array[j])
  {
   swap_temp = suffix_array[j-1];
   suffix_array[j-1] = suffix_array[j];
   suffix_array[j] = swap_temp;
   j--;
  }
  i++;
 }
 return suffix_array;
}

//inicializovanie pola positions, kazda hodnota urcuje zaciatocnu poziciu rotacie
//nasledne zotriedenie tohto pola, pricom sa porovnavaju rotacie retazca zacinajuce na tychto poziciach
unsigned int * mergesort_SA(unsigned int *suffix_array, char *s, unsigned int length){
  mergeSort(s,suffix_array,length);
  return suffix_array;
}

//procedure for constructing BWT of input string
unsigned char *create_bwt(unsigned int *suffix_array, unsigned char *s, unsigned int genome_length)
{
 unsigned int i;
 unsigned char *bwt = (unsigned char *) malloc (genome_length * sizeof(unsigned char));
 //suffix_array = insertion_sort_array(suffix_array,s,genome_length); 
 suffix_array = mergesort_SA(suffix_array,s,genome_length);
 for(i=0;i<genome_length;i++)
 {
  if (suffix_array[i]==0)
   bwt[i] = s[genome_length-1];
  else bwt[i] = s[suffix_array[i]-1];
 }
 return bwt;
}
