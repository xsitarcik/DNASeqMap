#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "compression.h"

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

char *move_to_front_decode (char *alphabet, unsigned int bitvector_length,unsigned char*bitvector)
{
 struct symbol_table *front = build_symbol_table(alphabet);
 struct symbol_table *current;
 int i, index;
 unsigned char shift, start;
 unsigned char bits_per_char = get_min_bits_per_char(alphabet);
 unsigned char code, next_code;
 int wrong_pos = wrong_pos = 8 - bits_per_char + 1;
 for (i=0;i<bitvector_length;i=i+bits_per_char)
 {
  index = i/8;
  start = i%8;
  code = bitvector[index];
  if (start<wrong_pos)
  {
   code = code >> (8 - start - bits_per_char);
  }
  else 
  {
   shift = start + bits_per_char - 8;
   code = code << shift;
   next_code = bitvector[index + 1];
   next_code = next_code >> (8 - shift);
   code = code | next_code;
  } 
  code = code & 7;
  current = decode(front,code);
  printf("%c",current->symbol);
  if (current->symbol != front->symbol)
   front = push_to_front(front,current);
 }
 return bitvector;
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

void move_to_front_encode (char *alphabet, char *s)
{
 unsigned char count;
 int i;
 struct symbol_table *front = build_symbol_table(alphabet);
 struct symbol_table *current = front;
 int string_length = strlen(s);
 //for each character in string
 for (i=0;i<string_length;i++)
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
}

unsigned char *arithmetic_encode (unsigned int *bitvector_length, char *s, char * alphabet, unsigned int string_length)
{
 unsigned char*bitvector = NULL;
 unsigned char remainder;
 unsigned char count;
 unsigned char bits_per_char = get_min_bits_per_char(alphabet);
 int i,j,k;
 int needed_bytes = bits_per_char*string_length/8 + 1;;
 int byte_index;
 int in_byte;
 int wrong_pos = 8 - bits_per_char + 1;
 printf("bits per char %d, string_length %d\n",bits_per_char, string_length);
 *bitvector_length = bits_per_char*string_length;
 bitvector = (unsigned char *)malloc(needed_bytes); 
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
  if (in_byte<wrong_pos)
  {
   bitvector[byte_index] = bitvector[byte_index] << bits_per_char;
   bitvector[byte_index] = bitvector[byte_index] + count;
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
  }
  j = j + bits_per_char;
 } 
 in_byte = j%8;
 byte_index = j/8;
 k = 8 - in_byte;
 bitvector[byte_index] = bitvector[byte_index] << k;
 
 return bitvector;
}

unsigned char get_index_in_alphabet(char *alphabet, char c)
{
 unsigned char i = 0;
 while(c!=alphabet[i])
  i++;
 return i;
}

void alphabet_encode (char *s, char * alphabet)
{
 unsigned char count;
 unsigned int string_length = strlen(s);
 unsigned int i;
 for (i=0;i<string_length;i++)
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
 int size = bitvector_length/8 + 1;
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
  printf(" ");
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


