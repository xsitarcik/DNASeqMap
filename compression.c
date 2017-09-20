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

unsigned char* move_to_front_encode (struct symbol_table *front, char *s)
{
 unsigned char*bitvector = NULL;
 unsigned char bits_per_char;
 int symbols_count = 1;
 int i = 2;
 unsigned int j = 1;
 int k;
 unsigned char count;
 unsigned char remainder;
 int string_length = strlen(s);
 int needed_bytes;
 int byte_index;
 int in_byte;
 int wrong_pos;
 struct symbol_table *current = front;
 struct symbol_table *swap;
 while (current->next!=NULL)
 {
  printf("inddex %d = %c\n",symbols_count,current->symbol);
  symbols_count++;
  current = current->next;
 }
 printf("inddex %d = %c\n",symbols_count,current->symbol);
 printf("pocet symbolov je %d\n",symbols_count);
 while (symbols_count>=i)
 {
  j++;
  i = i *2;
 }
 printf("we need %d bits per char\n", j);
 bits_per_char = j;
 wrong_pos = 8 - bits_per_char + 1;
 printf("totally we need %d bits\n",j*string_length);
 needed_bytes = j*string_length/8 + 1;
 printf("we'll be allocating %d bytes\n",needed_bytes);
 bitvector = (unsigned char *)malloc(needed_bytes); 
 if (bitvector==NULL)
 {
  printf("Error when allocating memory for bitvector in MTF\n");
  exit(1);
 }
 //for each character in string
 j = 0;
 for (i=0;i<string_length;i++)
 {
  count = 0;
  current = front;
  printf("na fronte je symbol %c\n",front->symbol);
  while (current->symbol != s[i])
  {
   current = current->next;
   count++;
  }
  printf("kod je %d, ukladam na j=%d\n",count,j);
  byte_index = j/8;
  in_byte = j%8;
  if (in_byte<wrong_pos)
  {
   bitvector[byte_index] = bitvector[byte_index] << bits_per_char;
   bitvector[byte_index] = bitvector[byte_index] + count;
   printf("na %d som supol %d\n",byte_index,count);
  }
  else
  {
   k = 8 - in_byte; //remaining bits in first byte
   remainder = count >> (bits_per_char - k);
   bitvector[byte_index] = bitvector[byte_index] << k;
   bitvector[byte_index] = bitvector[byte_index] + remainder;
   printf("na %d som supol %d\n",byte_index,remainder);
   k = get_set_bits(bits_per_char-k);
   printf("k je %d\n",k); 
   remainder = count & k;
   bitvector[byte_index + 1] = bitvector[byte_index + 1] + remainder;
   printf("na %d som supol %d\n",byte_index,remainder);
  }
  j = j + bits_per_char;
  if (count != 0)
  {
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
   front = current;
  } 
  while (current->next!=NULL)
   {
    printf("%c ",current->symbol); 
    current=current->next;
   }
   printf("%c\n",current->symbol);
  }
printf("\n");
 for (i=0;i<needed_bytes;i++)
  printf("%d ",bitvector[i]);
 return bitvector;
}

unsigned char get_set_bits(unsigned char bits)
{
 int i;
 int ret = 1;
 for (i=0;i<bits;i++)
  ret = ret * 2;
 return (ret - 1);
}

void print_bit_vector(unsigned char *bitvector)
{
 int size = strlen(bitvector);
 int i,j;
printf("\nsize je %d\n",size);
 for (i=0;i<size;i++)
 {
  for(j=0;j<8;j++)
  {
   printf("%i", bitvector[i] & 128 ? 1 : 0);
   bitvector[i] = bitvector[i] << 1;
  }
  printf(" ");
 }
}
