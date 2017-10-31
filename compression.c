#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "compression.h"


struct compressed_block *compress_FMIndex(unsigned int block_size, unsigned char flag_mtf, unsigned char flag_runs,unsigned char flag_huffman, 
	unsigned char*alphabet, unsigned char*bwt, unsigned int*length)
{

	// nutne rozbit bwt na bloky
 //vypocet poctu blokov
 unsigned char bits_per_char = get_min_bits_per_char(alphabet);
 printf("bits per %d\n",bits_per_char);
 unsigned int count_of_blocks = *length/block_size; //zvysok osetrit !
 unsigned int i,j,k,remainder;
 unsigned char *bitvector;
 unsigned int *temp_array;
 unsigned char*runs;
 unsigned int run_length;
 unsigned int temp_occ[strlen(alphabet)];
 struct compressed_block *array_of_blocks;
 array_of_blocks = (struct compressed_block*)malloc((count_of_blocks+1)*sizeof(struct compressed_block));
 j = 0;
 
 for (i = 0; i<count_of_blocks;i++)
 {
  array_of_blocks[i].occurences = (unsigned int*)malloc(strlen(alphabet)*sizeof(unsigned int));
  if (i!=0)
  {
    for (k=0;k<strlen(alphabet);k++)
     array_of_blocks[i].occurences[k] = array_of_blocks[i-1].occurences[k] + temp_array[k];
  }
  else 
  	for (k=0;k<strlen(alphabet);k++)
  	 array_of_blocks[i].occurences[k] = 0;
  temp_array = calc_occurences(&bwt[j],block_size,alphabet);
  
  array_of_blocks[i].block_size = block_size;
  printf("printing occ\n");
  for (k=0;k<strlen(alphabet);k++)
    printf("%d ",array_of_blocks[i].occurences[k]);
   printf("\n");

  if (flag_mtf)
  {
   printf("----Move to front encoding----\n");
   move_to_front_encode(alphabet,&bwt[j],block_size);
   printf("easz\n");
  }
  else
  {
   printf("----Alphabet Encoding----\n");
   alphabet_encode(&bwt[j],alphabet,block_size);
  }
  if (flag_runs)
  {
   run_length = block_size;
   printf("----Zero Runs Encoding----\n");
   printf("original length is %d\n",run_length);
   runs = run_length_encode(&bwt[j],&run_length);
   printf("new length after ZRE is %d\n",run_length);
  }
  
  if(flag_runs){
   array_of_blocks[i].bitvector = arithmetic_encode(runs,run_length,bits_per_char,&array_of_blocks[i].bitvector_length);
   free(runs);
  }
  else	
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
 	printf("handling remaining\n");
  count_of_blocks++;
  array_of_blocks[i].block_size = remainder;

  //temp_array = calc_occurences(&bwt[j],remainder,alphabet);
  array_of_blocks[i].occurences = (unsigned int*)malloc(strlen(alphabet)*sizeof(unsigned int));
  for (k=0;k<strlen(alphabet);k++)
   array_of_blocks[i].occurences[k] = array_of_blocks[i-1].occurences[k] + temp_array[k];
  printf("printing occ\n");
  for (k=0;k<strlen(alphabet);k++)
    printf("%d",array_of_blocks[i].occurences[k]);
   printf("\n");
  if (flag_mtf)
  {
   printf("----Move to front encoding----\n");
   move_to_front_encode(alphabet,&bwt[j],remainder);
  }
  else
  {
   printf("----Alphabet Encoding----\n");
   alphabet_encode(&bwt[j],alphabet,remainder);
  }
  if (flag_runs)
  {
   run_length = remainder;
   printf("----Zero Runs Encoding----\n");
   printf("original length is %d\n",run_length);
   runs = run_length_encode(&bwt[j],&run_length);
   printf("new length after ZRE is %d\n",run_length);
  }

  if (flag_runs)
  {
   array_of_blocks[i].bitvector = arithmetic_encode(runs,run_length,bits_per_char,&array_of_blocks[i].bitvector_length);
   free(runs);
  }
  else
   array_of_blocks[i].bitvector = arithmetic_encode(&bwt[j],remainder,bits_per_char,&array_of_blocks[i].bitvector_length);
  printf("bitvector length je %d\n",array_of_blocks[i].bitvector_length);
  print_bit_vector(array_of_blocks[i].bitvector, array_of_blocks[i].bitvector_length);
 }
 count_of_blocks--;
 printf("count of blocks je %d\n",count_of_blocks);
 *length = count_of_blocks;
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
 for (i=0;i<string_length;i++)
 {
  current = decode(front,s[i]);
  s[i]=current->symbol;
  if (current->symbol != front->symbol)
   front = push_to_front(front,current);
 }
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
 while (symbols_count>i)
 {
  j++;
  i = i *2;
 }
 return j;
}

void *move_to_front_encode (char *alphabet, char *s, unsigned int block_size)
{
 unsigned char count;
 int i;
 struct symbol_table *front = build_symbol_table(alphabet);
 struct symbol_table *current = front;
 printf("%s\n",s);
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

unsigned char *decompress_block(unsigned int bitvector_length, unsigned char*bitvector, unsigned int*runs, unsigned char flag_mtf, 
	unsigned char flag_runs, unsigned char flag_huffman, unsigned int block_size, unsigned char*alphabet)
{
 //printf("teraz dekodujem blok\n");
 unsigned int i;
 unsigned char *bwt = arithmetic_decode (bitvector,alphabet, block_size);
 //for (i=0;i<block_size;i++)
  //printf("%d",bwt[i]);
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
 /*printf("vypis\n");
 for(i=0;i<block_size;i++)
 	printf("%c",bwt[i]);*/
 return bwt;
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

unsigned char *arithmetic_decode (unsigned char *bitvector, char *alphabet, unsigned int string_length)
{
 unsigned char bits_per_char = get_min_bits_per_char(alphabet);
 //printf("velkost dekodu je %d\n",string_length);
 unsigned char*s = (unsigned char*)malloc(string_length);
 if (s==NULL)
 {
  printf("Error when allocating memory for bitvector\n");
  exit(1);
 }
 unsigned int i = 0;
 unsigned int bitposition = 0;
 while (i<string_length)
 {
  s[i] = decode_bits(bits_per_char,bitposition,bitvector);
  i++;
  bitposition = bitposition + bits_per_char;
 }
 //free(bitvector);
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

/*
RUN-LENGTH ENCODING
-when? = at least two same consecutive characters 
-what? = saves count of run-length of character (including those two)
-where? = in table "runs" of that block
*/
unsigned char* run_length_encode(unsigned char *s, unsigned int *string_length)
{
 unsigned int i = 0;
 unsigned int j = 0;
 unsigned int cmp = *string_length-1;
 unsigned char previous;
 unsigned char counter;
 unsigned char*encoded = (unsigned char *)malloc(*string_length*2);
 fflush(stdout);
 while(i<=cmp)
 {
  counter = 1;
  previous = s[i];
  encoded[j] = previous;
  j++;
  i++;  
  while (s[i]==previous && i<=cmp)
  {
   counter++;
   i++;
  }
  if (counter>1)
  {
  	printf("count je %d\n",counter);
   encoded[j] = previous;
   j++;
   encoded[j] = counter;
   j++;
  }
 }

 printf("encoded[j] %d\n",encoded[j-1]);
 printf("j = %d length = %d\n",j,*string_length); 
 if (j!=*string_length)
  encoded = (unsigned char *)realloc(encoded,j);
 printf("usetrilo sa %d bytov\n",*string_length-j);
 *string_length = j;
 printf("hodnota stringlength je %d\n",*string_length);
 //free(s);
 return encoded;
}

struct huffman_node *get_root(struct huffman_node **heap, unsigned int size)
{
    struct huffman_node* temp = heap[0];
    heap[0] = heap[size-1];
    heapify(heap, 0, size-1);
    return temp;
}

void heapify(struct huffman_node **heap, unsigned int index, unsigned int size)
{
 unsigned int left_child = index*2+1;
 unsigned int right_child = index*2+2;
 unsigned int min = index;
 struct huffman_node*swap;

 if (left_child<size && heap[left_child]->freq < heap[index]->freq)
  min = left_child;
 if (right_child<size && heap[right_child]->freq < heap[index]->freq)
  min = right_child;
 
 if (min != index)
 {
  swap = heap[min];
  heap[min] = heap[index];
  heap[index] = swap;
  heapify(heap,index,size);
 }
}

struct huffman_node *create_new_node(unsigned char c, unsigned int freq)
{
 struct huffman_node * current = (struct huffman_node *) malloc(sizeof(struct huffman_node));
 current->symbol = c;
 current->freq = freq;
 current->left = NULL;
 current->right = NULL;
 current->visited = 0;
 return current;
}

struct huffman_node **build_min_heap(unsigned char*alphabet, unsigned int*freq, unsigned int alphabet_length)
{
 int i;
 struct huffman_node ** heap = (struct huffman_node **) malloc(sizeof(struct huffman_node*)*alphabet_length);

 for (i=0;i<alphabet_length;i++)
 {
  heap[i] = create_new_node(alphabet[i],freq[i]);
 }
 for (i=0;i<alphabet_length;i++)
 {
 	printf("%d=%d\n",heap[i]->freq,heap[i]->symbol);
 }
 printf("\n");
 for (i=alphabet_length/2;i>=0;i--)
 {
 	heapify(heap,i,alphabet_length);
 }
for (i=0;i<alphabet_length;i++)
 {
 	printf("%d=%d\n",heap[i]->freq,heap[i]->symbol);
 }
 printf("\n");

 return heap;
}

struct huffman_node *build_huffman_tree(unsigned char*alphabet,unsigned int*freq)
{
 unsigned int alphabet_length = strlen(alphabet);
 struct huffman_node ** heap = build_min_heap(alphabet,freq,alphabet_length);
 struct huffman_node *left_node;
 struct huffman_node *right_node;
 struct huffman_node *root_node;
 struct huffman_node *current;
 unsigned int i = 0;
 unsigned int j = 0;
 int size = alphabet_length;
 unsigned char code[alphabet_length];
 unsigned int huffman_tree_size = size+size-1;
 while (size>1)
 {
    printf("size je %d\n",size);
 //extract 2 min frequencies from heap
 left_node = get_root(heap,size);
 printf("extrahovali sme %d %c\n",left_node->freq,left_node->symbol);
 size--;
 if (size!=0)
 {
  right_node = get_root(heap,size);
  size--;
 }
 else
 {
  right_node = heap[0];
  size--;
 }
printf("extrahovali sme %d %c\n",right_node->freq,right_node->symbol);
 //add extracted frequencies and its sum to huffman tree
 //add to heap sum of extracted roots of minheap
 root_node = create_new_node('a',left_node->freq+right_node->freq);
 root_node->left = left_node;
 root_node->right = right_node;
 heap[size] = root_node;
 size++;
 }
 printf("ukonecne\n");

 
 while (i<alphabet_length)
 {
  current = root_node;
  j = 0;
  while(1)
  {
   while (current->left != NULL && current->left->visited==0)
   {
    current = current->left;
    code[j] = '0';
    j++;
   }
   while (current->right != NULL && current->right->visited==0)
    {
    current = current->right;
    code[j] = '1';
    j++;
   }
   //ak je list
   if (current->left==NULL && current->right==NULL)
   {
    current->visited = 1;
    break;
   }
   else if (current->left->visited==1 && current->right->visited==1)
   {
    current->visited = 1;
    break;
   }
  }
  //ak je list
   if (current->left==NULL && current->right==NULL)
   {
    code[j]='\0';
    printf("kod je %s, %d %c\n",code,current->freq,current->symbol);
    i++;
   }
 }
 printf("\n");
 
}


