#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "compression.h"


struct compressed_block *compress_FMIndex(unsigned int block_size, unsigned char flag_mtf, unsigned char flag_zero_runs,unsigned char flag_huffman, 
	unsigned char*alphabet, unsigned char*bwt, unsigned int*length)
{

	// nutne rozbit bwt na bloky
 //vypocet poctu blokov
 unsigned char bits_per_char = get_min_bits_per_char(alphabet);;
 unsigned int count_of_blocks = *length/block_size; //zvysok osetrit !
 unsigned int i,j,k, remainder;
 unsigned char *bitvector;
 struct compressed_block *array_of_blocks;
 array_of_blocks = (struct compressed_block*)malloc((count_of_blocks+1)*sizeof(struct compressed_block));
 j = 0;
 
 for (i = 0; i<count_of_blocks;i++)
 {
  array_of_blocks[i].occurences = calc_occurences(&bwt[j],block_size,alphabet);
  printf("printing occ\n");
  for (k=0;k<strlen(alphabet);k++)
    printf("%d",array_of_blocks[i].occurences[k]);
   printf("\n");

  if (flag_mtf)
  {
   printf("----Move to front encoding----\n");
   array_of_blocks[i].front = move_to_front_encode(alphabet,&bwt[j],block_size);
  }
  else
  {
   printf("----Alphabet Encoding----\n");
   alphabet_encode(&bwt[j],alphabet,block_size);
  }
  if (flag_zero_runs)
  {
  	//TO TREBA PREROBIT
  	//alloc runs
  printf("----Zero Runs Encoding----\n");
  printf("original length is %d\n",*length);
  bwt = zero_runs_encode(bwt,length);
  printf("new length is %d\n",*length);
  //array_of_blocks[i].runs = 0;
  }
  
  array_of_blocks[i].bitvector = arithmetic_encode(&bwt[j],block_size,bits_per_char,&array_of_blocks[i].bitvector_length);
  printf("bitvector length je %d\n",array_of_blocks[i].bitvector_length);
  print_bit_vector(array_of_blocks[i].bitvector, array_of_blocks[i].bitvector_length);
   
  j = j + block_size;
 }

 //HANDLING REMAINING BWT
 printf("zostava %d\n",*length - j);
 remainder = *length - j;
 if (remainder>0)
 {
  array_of_blocks[i].occurences = calc_occurences(&bwt[j],remainder,alphabet);
  printf("printing occ\n");
  for (k=0;k<strlen(alphabet);k++)
    printf("%d",array_of_blocks[i].occurences[k]);
   printf("\n");
  if (flag_mtf)
  {
   printf("----Move to front encoding----\n");
   array_of_blocks[i].front = move_to_front_encode(alphabet,&bwt[j],remainder);
  }
  else
  {
   printf("----Alphabet Encoding----\n");
   array_of_blocks[i].front = NULL;
   alphabet_encode(&bwt[j],alphabet,remainder);
  }
  if (flag_zero_runs)
  {
  	//TO TREBA PREROBIT
  	//alloc runs
  printf("----Zero Runs Encoding----\n");
  printf("original length is %d\n",*length);
  bwt = zero_runs_encode(bwt,length);
  printf("new length is %d\n",*length);
  //array_of_blocks[i].runs = 0;
  }

  bits_per_char = get_min_bits_per_char(alphabet);
  array_of_blocks[i].bitvector = arithmetic_encode(&bwt[j],remainder,bits_per_char,&array_of_blocks[i].bitvector_length);
  printf("bitvector length je %d\n",array_of_blocks[i].bitvector_length);
  print_bit_vector(array_of_blocks[i].bitvector, array_of_blocks[i].bitvector_length);
 }


 free(bwt);

  

 return array_of_blocks;
}

unsigned int *calc_occurences(char *s, int string_length, char *alphabet)
{
 unsigned int alphabet_size = strlen(alphabet);
 unsigned int *occurence_table = (unsigned int *)calloc(alphabet_size,sizeof(unsigned int));
 unsigned int i,j;
 if (occurence_table==NULL)
  {
   printf("Error. No memory when allocating occurence table\n");
   exit(1);
  }

 for (i=0;i<string_length;i++)
 {
  for (j=0;j<alphabet_size;j++)
   if (s[i]==alphabet[j])
    occurence_table[j]++;     
 } 
 return occurence_table;
}

struct symbol_table *build_symbol_table(char *alphabet)
{
 int alphabet_size = strlen(alphabet);
 unsigned int i = 0;
 struct symbol_table *front = NULL;
 struct symbol_table *prev = (struct symbol_table*)malloc(sizeof(struct symbol_table));
 struct symbol_table *next = NULL;
 prev->symbol = alphabet[i];
 prev->previous = NULL;
 if (prev==NULL)
  {
   printf("Error when allocating memory for symbol table\n");
   exit(1);
  }
 front = prev; //for returning
 for (i=1;i<alphabet_size;i++){
  next = (struct symbol_table*)malloc(sizeof(struct symbol_table));
  if (next==NULL)
  {
   printf("Error when allocating memory for symbol table\n");
   exit(1);
  }
  next->symbol = alphabet[i];
  next->previous = prev;
  prev->next = next;
  prev = next;
 }
 next->next = NULL;
 return front;
}

struct symbol_table *push_to_front(struct symbol_table *front, struct symbol_table *current)
{
 struct symbol_table *swap;
 if (current->previous!=NULL)
 {
  swap = current->previous;
  swap->next = current->next;
 }
 if (current->next!=NULL)
 {
  swap = current->next;
  swap->previous = current->previous;
 }
 current->previous = NULL;
 current->next = front;
 front->previous = current;
 return current;
}

void *move_to_front_decode (char *alphabet, unsigned int string_length,unsigned char*s)
{
 struct symbol_table *front = build_symbol_table(alphabet);
 struct symbol_table *current;
 unsigned int i;
 printf("string_length = %d\n",string_length);
 for (i=0;i<string_length;i++)
 {
  printf("%d,s[i] e %d = ",i,s[i]);
  current = decode(front,s[i]);
  s[i]=current->symbol;
  printf("%c %c\n",current->symbol, front->symbol);
  if (current->symbol != front->symbol)
   front = push_to_front(front,current);
  printf("gr\n");
 }
 printf("wut\n");
 printf("wut\n");
 printf("wut\n");
 printf("wut\n");
 printf("wut\n");
}

struct symbol_table *decode(struct symbol_table *front, unsigned char code)
{
 unsigned char i = 0;
 struct symbol_table *current = front;
 while (i!=code)
 {
  current = current->next;
  i++;
 }
 return current;
}

unsigned char get_min_bits_per_char(char *alphabet)
{
unsigned int symbols_count = strlen(alphabet);
unsigned int i = 2;
unsigned int j = 1;
 while (symbols_count>=i)
 {
  j++;
  i = i *2;
 }
 return j;
}

struct symbol_table *move_to_front_encode (char *alphabet, char *s, unsigned int block_size)
{
 unsigned char count;
 int i;
 struct symbol_table *front = build_symbol_table(alphabet);
 struct symbol_table *current = front;
 //for each character in string
 for (i=0;i<block_size;i++)
 {
  count = 0;
  current = front;
  while (current->symbol != s[i])
  {
   current = current->next;
   count++;
  }
  s[i] = count;
  printf("%d",s[i]);
  if (count != 0)
   front = push_to_front(front,current);
 }
 printf("\n");
 return front;
}

unsigned char *arithmetic_encode (char *s, unsigned int string_length, unsigned char bits_per_char, unsigned int *bitvector_length)
{
 unsigned char*bitvector = NULL;
 unsigned char remainder;
 unsigned char count;
 int i,j,k;
 int needed_bytes = (string_length*bits_per_char+7)/8;
 //printf("needed bytes is %d\n",needed_bytes);
 int byte_index;
 int in_byte;
 int wrong_pos = 8 - bits_per_char + 1;
 //printf("bits per char %d, string_length %d\n",bits_per_char, string_length);
 unsigned int size = bits_per_char*string_length;
 bitvector = (unsigned char *)calloc(needed_bytes,1); 
 if (bitvector==NULL)
 {
  printf("Error when allocating memory for bitvector\n");
  exit(1);
 }
 //for each character in string
 j = 0;
 for (i=0;i<string_length;i++)
 {
  count = s[i];
  byte_index = j/8;
  in_byte = j%8;
  //printf("znak je %d, byte_index %d\n",s[i],byte_index);
  if (in_byte<wrong_pos)
  {
   bitvector[byte_index] = bitvector[byte_index] << bits_per_char;
   bitvector[byte_index] = bitvector[byte_index] + count;
   //printf("stav je %d\n",bitvector[byte_index]);
  }
  else 
  {
   k = 8 - in_byte; //remaining bits in first byte
   remainder = count >> (bits_per_char - k);
   bitvector[byte_index] = bitvector[byte_index] << k;
   bitvector[byte_index] = bitvector[byte_index] + remainder;
   k = get_set_bits(bits_per_char-k);
   remainder = count & k;
   bitvector[byte_index + 1] = bitvector[byte_index + 1] + remainder;
   //printf("stav je %d a dalsi %d\n",bitvector[byte_index],bitvector[byte_index+1]);
  }
  j = j + bits_per_char;
 } 
 in_byte = j%8;
 byte_index = j/8;
 k = 8 - in_byte;
 bitvector[byte_index] = bitvector[byte_index] << k;
 *bitvector_length = size;
 return bitvector;
}

unsigned char *decompress_block(unsigned int bitvector_length, struct symbol_table *front, unsigned char*bitvector, unsigned int*runs, unsigned char flag_mtf, 
	unsigned char flag_zero_runs, unsigned char flag_huffman, unsigned int block_size, unsigned char*alphabet)
{
 printf("teraz deokduujem\n");
 unsigned int i;
 unsigned char *bwt = arithmetic_decode (bitvector,alphabet,bitvector_length, block_size);
 for (i=0;i<block_size;i++)
  printf("%d",bwt[i]);
 if (flag_mtf)
 {
  printf("------MTF decoding----\n");
  move_to_front_decode(alphabet, block_size, bwt);
 }
 else
 {
  for(i=0;i<block_size;i++)
  {
   bwt[i] = alphabet[bwt[i]];
  }
 }
 for(i=0;i<block_size;i++)
 	printf("%s\n",bwt);
} 

unsigned char decode_bits(unsigned char bits_per_char,unsigned int bitposition, unsigned char *bitvector)
{
 unsigned char remainder = bitposition%8;
 unsigned char maxpos = 8 - bits_per_char;
 unsigned int index = bitposition/8;
 unsigned char result = bitvector[index];
 unsigned char temp;
 //check if code is contained in 2 bytes
 if (remainder>maxpos)
 {
  unsigned char bits_in_first_byte = 8-remainder;
  unsigned char bits_in_second_byte = bits_per_char - bits_in_first_byte;
  result = result & get_set_bits(bits_in_first_byte);
  result = result << (bits_in_second_byte);
  temp = bitvector[index+1] >> (8 - bits_in_second_byte);
  temp = temp & get_set_bits(bits_in_second_byte);
  result = result + temp;
 }
 else
 {
  result = result >> (maxpos - remainder);
  result = result & get_set_bits(bits_per_char);
 }
 return result;
}

unsigned char *arithmetic_decode (unsigned char *bitvector, char *alphabet, unsigned int bitvector_length, unsigned int string_length)
{
 unsigned char bits_per_char = get_min_bits_per_char(alphabet);
 printf("velkost dekodu je %d\n",string_length);
 unsigned char*s = (unsigned char*)malloc(string_length);
 if (s==NULL)
 {
  printf("Error when allocating memory for bitvector\n");
  exit(1);
 }
 unsigned int i = 0;
 unsigned int bitposition = 0;
 while (bitposition<bitvector_length)
 {
  s[i] = decode_bits(bits_per_char,bitposition,bitvector);
  i++;
  bitposition = bitposition + bits_per_char;
 }
 free(bitvector);
 return s;
}

unsigned char get_index_in_alphabet(char *alphabet, char c)
{
 unsigned char i = 0;
 while(c!=alphabet[i])
  i++;
 return i;
}

void alphabet_encode (char *s, char * alphabet,unsigned int block_size)
{
 unsigned char count;
 unsigned int i;
 for (i=0;i<block_size;i++)
 {
  count = get_index_in_alphabet(alphabet,s[i]);
  s[i] = count;
  printf("%d",s[i]);
 }
 printf("\n");
}

unsigned char get_set_bits(unsigned char bits)
{
 int i;
 int ret = 1;
 for (i=0;i<bits;i++)
  ret = ret * 2;
 return (ret - 1);
}

void print_bit_vector(unsigned char *bitvector, unsigned int bitvector_length)
{
 int i,j;
 int size = (bitvector_length+7)/8;
 unsigned char code;
 printf("Actual length is: %d bits\n",bitvector_length);
 for (i=0;i<size;i++)
 {
  code = bitvector[i];
  for(j=0;j<8;j++)
  {
   printf("%i", code & 128 ? 1 : 0);
   code = code << 1;
  }
 }
 printf("\n");
}

unsigned char* zero_runs_encode(unsigned char *s, unsigned int *string_length)
{
 unsigned int i = 0;
 unsigned int j = 0;
 unsigned int cmp = *string_length-1;
 unsigned char *encoded = (unsigned char*) malloc(cmp*2*sizeof(char));
 unsigned char previous = s[0];
 unsigned char counter = 0;
 printf("hodnota stringlenght je %d\n",*string_length);
 while(i<=cmp)
 {
  counter = 0;  
  while (s[i]==0 && i!=cmp)
  {
   counter++;
   i++;
  }
  if (counter!=0)
  {
   encoded[j] = 0;
   j++;
   encoded[j] = counter;
   j++;
  }
  else 
  {
   encoded[j] = s[i];
   j++;
   i++;
  }
 }
 printf("encoded[j] %d\n",encoded[j-1]);
 printf("finne %d %d\n",j,*string_length); 
 if (j!=*string_length)
  encoded = (unsigned char *)realloc(encoded,j);
 printf("usetrilo sa %d bytov\n",*string_length-j);
 *string_length = j;
 printf("hodnota stringlenght je %d\n",*string_length);
 free(s);
 return encoded;
}


