#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "compression.h"

//procedure dividing BWT string into array of blocks and compressing each block
struct compressed_block *compress_FMIndex(unsigned int block_size, unsigned char flag_mtf, unsigned char flag_runs,unsigned char flag_huffman, 
	unsigned char*alphabet, unsigned char*bwt, unsigned int*length)
{
 unsigned char bits_per_char = get_min_bits_per_char(alphabet);
 unsigned int count_of_blocks = *length/block_size; //get count of needed blocks
 unsigned int i,j,k,remainder;
 unsigned char *bitvector;
 unsigned int *temp_array;
 unsigned char*runs;
 unsigned int run_length;
 unsigned char*temp_alphabet;
 unsigned char*new_alphabet;
 unsigned int *frequencies;
 struct compressed_block *array_of_blocks;
 unsigned char new_alphabet_length;
 unsigned char alphabet_index;
 unsigned int temp_alphabet_length;
 array_of_blocks = (struct compressed_block*)malloc((count_of_blocks+1)*sizeof(struct compressed_block));
 j = 0;

 //process each block
 for (i = 0; i<count_of_blocks;i++)
 {
  
  //create table of occurences for block
  array_of_blocks[i].occurences = (unsigned int*)malloc(strlen(alphabet)*sizeof(unsigned int));
  
  //table of occurences for first block should be zeroes, else count of previous block + occurences in previous block
  if (i!=0)
  {
    for (k=0;k<strlen(alphabet);k++)
     array_of_blocks[i].occurences[k] = array_of_blocks[i-1].occurences[k] + temp_array[k];
  }
  else 
  	for (k=0;k<strlen(alphabet);k++)
  	 array_of_blocks[i].occurences[k] = 0;

  //calculate occurences for this block
  temp_array = calc_occurences(&bwt[j],block_size,alphabet);
  
  //set block_size for block, so we wont access unallocated memory
  array_of_blocks[i].block_size = block_size;

  //handle move-to-front transform if needed
  if (flag_mtf)
  {
   move_to_front_encode(alphabet,&bwt[j],block_size);
  }
  else
  {
    //else we need to transform from alphabetical to numerical A,C,G,T to 0,1,2,3 for example
   alphabet_encode(&bwt[j],alphabet,block_size);
  }

  //handle run-length encoding if needed
  if (flag_runs)
  {

   //prepare auxiliary variable for getting size of RLE encoded string
   run_length = block_size;

   //prepare auxiliary string of new alphabet added when RLE encoding
   temp_alphabet = (unsigned char*)malloc(block_size/2);

   //encode block of BWT string, send auxiliary variables as well
   runs = run_length_encode(&bwt[j],&run_length,temp_alphabet,&temp_alphabet_length,flag_runs);

   //realloc auxiliary temp alphabet with real size of alphabet
   temp_alphabet = (unsigned char*)realloc(temp_alphabet,temp_alphabet_length);
   new_alphabet = order_new_alphabet(temp_alphabet,alphabet,&new_alphabet_length,temp_alphabet_length);

   for(k=0;k<strlen(new_alphabet);k++)
    printf("pocet: %d\n",new_alphabet[k]);
   printf("new length after ZRE is %d\n",run_length);
  }

  //handle huffman encoding if needed
  if(flag_huffman)
  {
   
   //check if input is RLE encoded
   if(flag_runs)
   {
    //from constructed new alphabet get frequencies
    /*printf("nova abeceda je dlzky %d a:\n",new_alphabet_length);
    for (k=0;k<new_alphabet_length;k++)
     printf("%d",new_alphabet[k]);*/

    //count frequencies of characters in alphabet
    frequencies = (unsigned int*)calloc(new_alphabet_length,sizeof(unsigned int));
    for (k=0;k<run_length;k++)
    {
     if(runs[k]>3)
        printf("hladam index %d ajaj %d\n",k,runs[k]);
        //printf("wut %d %d\n",k,runs[k]);
     alphabet_index = get_index_in_alphabet(new_alphabet,runs[k]);
     frequencies[alphabet_index]++;

    }

    //build huffman tree and save its root to block
    array_of_blocks[i].huffman_tree = build_huffman_tree(new_alphabet, frequencies,new_alphabet_length);

    //huffman encode block of BWT string with help of huffman tree
    array_of_blocks[i].bitvector= pack_huffman_to_bitvector(array_of_blocks[i].huffman_tree,new_alphabet,new_alphabet_length,runs,run_length,&array_of_blocks[i].bitvector_length,1);
    free(runs);
   }
   else
   {
    new_alphabet = (unsigned char*)malloc(strlen(alphabet));
    for (k=0;k<strlen(alphabet);k++)
     new_alphabet[k]=k;
    frequencies = (unsigned int*)calloc(strlen(alphabet),sizeof(unsigned int));
    for (k=0;k<block_size;k++)
    {
     frequencies[bwt[j+k]]++;
    }
    array_of_blocks[i].huffman_tree = build_huffman_tree(new_alphabet, frequencies,strlen(alphabet));
    array_of_blocks[i].bitvector= pack_huffman_to_bitvector(array_of_blocks[i].huffman_tree,new_alphabet,strlen(alphabet),&bwt[j],block_size,&array_of_blocks[i].bitvector_length,0);
   }
   printf("bitvector length je %d\n",array_of_blocks[i].bitvector_length);
   free(frequencies);
   //free(new_alphabet);
  }
  else if(flag_runs)
  {
   array_of_blocks[i].bitvector = bit_pack(runs,run_length,bits_per_char,&array_of_blocks[i].bitvector_length);
   free(runs);
  }
  else	
   array_of_blocks[i].bitvector = bit_pack(&bwt[j],block_size,bits_per_char,&array_of_blocks[i].bitvector_length);
  
  printf("bitvector length je %d j je %d\n",array_of_blocks[i].bitvector_length,j);
  //print_bit_vector(array_of_blocks[i].bitvector, array_of_blocks[i].bitvector_length);
  j = j + block_size;
 }

 //HANDLING LAST REMAINING BLOCK
 //get size of remaining block
 printf("handling remaining\n");
 remainder = *length - j;
 if (remainder>0)
 {
  count_of_blocks++;
  array_of_blocks[i].block_size = remainder;

  //calculate occurences
  array_of_blocks[i].occurences = (unsigned int*)malloc(strlen(alphabet)*sizeof(unsigned int));
  for (k=0;k<strlen(alphabet);k++)
   array_of_blocks[i].occurences[k] = array_of_blocks[i-1].occurences[k] + temp_array[k];
  
  printf("wtf remainder je %d j je %d\n",remainder,j);
  fflush(stdout);
  //handle MTF if needed
  if (flag_mtf)
  {
   move_to_front_encode(alphabet,&bwt[j],remainder);
  }
  else
  {
   alphabet_encode(&bwt[j],alphabet,remainder);
  }

  printf("wtf\n");
  fflush(stdout);

  //handle RLE if needed
  if (flag_runs)
  {
   run_length = remainder;
   //printf("----Zero Runs Encoding----\n");
   //printf("original length is %d\n",run_length);
   temp_alphabet = (unsigned char*)malloc(remainder/2);
   runs = run_length_encode(&bwt[j],&run_length,temp_alphabet,&temp_alphabet_length,flag_runs);
   /*for (k=0;k<run_length;k++)
    printf("%d",runs[k]);
   printf("\n");*/
   new_alphabet = order_new_alphabet(temp_alphabet,alphabet,&new_alphabet_length,temp_alphabet_length);
   /*for(k=0;k<strlen(new_alphabet);k++)
    printf("pocet: %d\n",new_alphabet[k]);
   printf("new length after ZRE is %d\n",run_length);*/
  }

  if(flag_huffman)
  {
   //printf("----Huffman Encoding----\n");
   if (flag_runs)
   {
    frequencies = (unsigned int*)calloc(new_alphabet_length,sizeof(unsigned int));
    for (k=0;k<run_length;k++)
    {
     alphabet_index = get_index_in_alphabet(new_alphabet,runs[k]);
     frequencies[alphabet_index]++;
    }
    array_of_blocks[i].huffman_tree = build_huffman_tree(new_alphabet, frequencies, new_alphabet_length);
    array_of_blocks[i].bitvector= pack_huffman_to_bitvector(array_of_blocks[i].huffman_tree,new_alphabet,new_alphabet_length, runs,run_length,&array_of_blocks[i].bitvector_length,1);
    //free(runs);
   }
   else
   {
    new_alphabet = (unsigned char*)malloc(strlen(alphabet));
    for (k=0;k<strlen(alphabet);k++)
     new_alphabet[k]=k;
    frequencies = (unsigned int*)calloc(strlen(alphabet),sizeof(unsigned int));
    for (k=0;k<remainder;k++)
    {
     frequencies[bwt[j+k]]++;
    }
    array_of_blocks[i].huffman_tree = build_huffman_tree(new_alphabet, frequencies,strlen(alphabet));
    array_of_blocks[i].bitvector= pack_huffman_to_bitvector(array_of_blocks[i].huffman_tree,new_alphabet,strlen(alphabet),&bwt[j],remainder,&array_of_blocks[i].bitvector_length,0);
   }
   printf("HF bitvector length je %d\n",array_of_blocks[i].bitvector_length);
   //free(frequencies);
   //free(new_alphabet);
  }
  else if(flag_runs)
  {
   array_of_blocks[i].bitvector = bit_pack(runs,run_length,bits_per_char,&array_of_blocks[i].bitvector_length);
   //free(runs);
  }
  else
   array_of_blocks[i].bitvector = bit_pack(&bwt[j],remainder,bits_per_char,&array_of_blocks[i].bitvector_length);
  printf("HFbitvector length je %d\n",array_of_blocks[i].bitvector_length);
  //print_bit_vector(array_of_blocks[i].bitvector, array_of_blocks[i].bitvector_length);
 }
 count_of_blocks--;
 printf("count of blocks je %d\n",count_of_blocks);
 *length = count_of_blocks;
 //free(bwt);

  

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
  /*printf("%d",s[i]);
  fflush(stdout);*/
  if (count != 0)
   front = push_to_front(front,current);
 }
 //printf("\n");
 return front;
}

unsigned char *bit_pack (char *s, unsigned int string_length, unsigned char bits_per_char, unsigned int *bitvector_length)
{
 unsigned char*bitvector = NULL;
 unsigned char remainder;
 unsigned char count;
 int i,j,k;
 int needed_bytes = (string_length*bits_per_char+7)/8;
 printf("needed bytes is %d\n",needed_bytes);
 int byte_index;
 int in_byte;
 int wrong_pos = 8 - bits_per_char + 1;
 printf("bits per char %d, string_length %d\n",bits_per_char, string_length);
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

unsigned char *decompress_block(unsigned int bitvector_length, unsigned char*bitvector, unsigned char flag_mtf, 
	unsigned char flag_runs, unsigned char flag_huffman, unsigned int block_size, unsigned char*alphabet, struct huffman_node*huffman_tree)
{
 //printf("teraz dekodujem blok\n");
 unsigned int i;
 unsigned char *bwt;
 if (flag_huffman)
 {
  if (flag_runs)
   bwt = huffman_decode_with_RLE(bitvector, huffman_tree, block_size,flag_runs);
  else
   bwt = huffman_decode_without_RLE(bitvector, huffman_tree, block_size);
 }
 else
 {
  if (flag_runs)
   bwt = arithmetic_decode_with_RLE(bitvector,alphabet, block_size);
  else
   bwt = arithmetic_decode_without_RLE(bitvector,alphabet, block_size);
 }

 if (flag_mtf)
 {
  //printf("\n------MTF decoding----\n");
  move_to_front_decode(alphabet, block_size, bwt);
  //printf("\n------MTF decoding succesfull----\n");
 }
 
 else
 {
  for(i=0;i<block_size;i++)
  {
   bwt[i] = alphabet[bwt[i]];
  }
 }
 /*
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

unsigned char *arithmetic_decode_with_RLE (unsigned char *bitvector, char *alphabet, unsigned int string_length)
{
 unsigned char bits_per_char = get_min_bits_per_char(alphabet);
 unsigned char repeat;
 //printf("velkost dekodu je %d\n",string_length);
 unsigned char*s = (unsigned char*)calloc(string_length,1);
 if (s==NULL)
 {
  printf("Error when allocating memory for bitvector\n");
  exit(1);
 }
 printf("HUE\n");
 unsigned int i = 0;
 unsigned char j;
 unsigned int bitposition = 0;
 unsigned char is_run = 0;
 
 //GET FIRST SYMBOL
 s[i] = decode_bits(bits_per_char,bitposition,bitvector);
 i++;
 bitposition = bitposition + bits_per_char;
 
 while (i<string_length)
 {
  s[i] = decode_bits(bits_per_char,bitposition,bitvector);
  bitposition = bitposition + bits_per_char;
  if (s[i]==s[i-1])
  {
   j = 0;
   repeat = decode_bits(bits_per_char,bitposition,bitvector);
   bitposition = bitposition + bits_per_char;
   while (i<string_length && j<repeat)
   {
    i++;
    s[i]=s[i-1];
    j++;
   }
  }
   i++;
 }

 //free(bitvector);
 return s;
}

unsigned char *arithmetic_decode_without_RLE (unsigned char *bitvector, char *alphabet, unsigned int string_length)
{
 unsigned char bits_per_char = get_min_bits_per_char(alphabet);
 //printf("velkost dekodu je %d\n",string_length);
 unsigned char*s = (unsigned char*)malloc(string_length);
 if (s==NULL)
 {
  printf("Error when allocating memory for bitvector\n");
  exit(1);
 }
 //printf("-----arithmetic decoding without RLE-----\n");
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

unsigned char* huffman_decode_without_RLE(unsigned char*bitvector, struct huffman_node*huffman_tree,unsigned int string_length)
{
 unsigned char*s = (unsigned char*)malloc(string_length);
 struct huffman_node *current = huffman_tree;
 unsigned int i=0;
 unsigned int bit_position = 8;
 unsigned int byte_index = 0;
 if (s==NULL)
 {
  printf("Error when allocating memory for bitvector\n");
  exit(1);
 }
 
 //printf("-----Huffman decoding without RLE-----\n");

 while (i<string_length)
 {
  //start from root of huffman tree
  current = huffman_tree;

  //until leaf is reached traverse tree according to code
  while(current->left!=NULL && current->right!=NULL)
  {
   if (CHECK_BIT(bitvector[byte_index],--bit_position))
    current = current->right;
   else 
    current = current->left;
   if (bit_position==0)
   {
    bit_position = 8;
    byte_index++;
   }
  }
  //save decoded code
  s[i] = current->symbol;
  //printf("%d",s[i]);
  i++;
 }
 return s;
}

unsigned char* huffman_decode_with_RLE(unsigned char*bitvector, struct huffman_node*huffman_tree,unsigned int string_length, unsigned char flag_runs)
{
 unsigned char*s = (unsigned char*)calloc(string_length+1,1);
 struct huffman_node *current = huffman_tree;
 unsigned int i=0;
 unsigned char j;
 unsigned int bit_position = 8;
 unsigned int byte_index = 0;
 unsigned char is_run = 0;
 if (s==NULL)
 {
  printf("Error when allocating memory for bitvector\n");
  exit(1);
 }
 
 //printf("-----Huffman decoding with RLE %d-----\n",string_length);
 
 //GET FIRST CHAR
 //until leaf is reached traverse tree according to code
  while(current->left!=NULL && current->right!=NULL)
  {
   if (CHECK_BIT(bitvector[byte_index],--bit_position))
    current = current->right;
   else 
    current = current->left;
   if (bit_position==0)
   {
    bit_position = 8;
    byte_index++;
   }
  }
  //save decoded code
  s[i] = current->symbol;
  //printf("%d",s[i]);
  i++;
 
 while (i<string_length)
 {
  //start from root of huffman tree
  current = huffman_tree;

  //until leaf is reached traverse tree according to code
  while(current->left!=NULL && current->right!=NULL)
  {
   if (CHECK_BIT(bitvector[byte_index],--bit_position))
    current = current->right;
   else 
    current = current->left;
   if (bit_position==0)
   {
    bit_position = 8;
    byte_index++;
   }
  }
  if (is_run>=flag_runs)
  {
   j = 0;

   //printf("bude %d opakovani, i je %d\n",current->symbol,i);
   while (j<current->symbol && i<string_length)
   {
    s[i]=s[i-1];
    //printf("%d",s[i]);
    i++;
    j++;
   }
   if (j==current->symbol)
    is_run = 0;
   //printf("j je %d opakovani, i je %d\n",j,i);
  }
  else
  {
   //save decoded code
   s[i] = current->symbol;
   //printf("%d",s[i]);
   if (s[i]==s[i-1])
    is_run++;
   else
    is_run = 0;
   i++;
  }
 }
 //printf("i je %d\n",i);
 s[i]='\0';
 //printf("%s",s);
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
  //printf("%d",s[i]);
 }
 //printf("\n");
}

unsigned char get_set_bits(unsigned char bits)
{
 int ret = 1;
 for (unsigned int i=0;i<bits;i++)
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
10122322211010030333120310000030012231323032100220203131201010201212202212131122013220001013233212332003101313330231031122102010
10122032211100100030331120310033000122031323032100022002031312010102012122002201213110220013220001101323302123302000310131331023103110220102010

10122322211010030333120310000030012231323032100220203131201010201212202212131122013220001013233212332003101313330231031122102301


02301231222233101200310113120012121231213301303111033212101301300302201013021021313130012023132110033102131302231331
023012312223301012000310110312000121212312133001303111033021210130130003022001013021021313130001202313211000033010213130220313301
02301231222233101200310113120012121231213301303111033212101301300302201013021021313130012023132110033102131302231332
*/

/*
RUN-LENGTH ENCODING
-when? = at least two same consecutive characters 
-what? = saves count of run-length of character (including those two)
-where? = in table "runs" of that block
*/
/*
unsigned char* run_length_encode(unsigned char *s, unsigned int *string_length, unsigned char* new_alphabet, unsigned int *new_alphabet_length)
{
 unsigned int i = 0;
 unsigned int j = 0;
 unsigned int k = 0;
 unsigned int cmp = *string_length-1;
 unsigned char previous;
 unsigned char counter;
 unsigned char*encoded = (unsigned char *)calloc(*string_length*2,1);
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
   //printf("%d,",counter);
   encoded[j] = previous;
   j++;
   encoded[j] = counter-2;
   j++;
   new_alphabet[k] = counter-2;
   k++;
  }
 }

 printf("j = %d length = %d\n",j,*string_length); 
 if (j!=*string_length)
  encoded = (unsigned char *)realloc(encoded,j);
 printf("usetrilo sa %d bytov\n",*string_length-j);
 *string_length = j;
 printf("hodnota string_length je %d\n",*string_length);
 *new_alphabet_length = k;

 //free(s);
 return encoded;
}*/


unsigned char* run_length_encode(unsigned char *s, unsigned int *string_length, unsigned char* new_alphabet, unsigned int *new_alphabet_length, unsigned char flag_runs)
{
 unsigned int i = 0;
 unsigned int j = 0;
 unsigned int k = 0;
 unsigned char c = 0;
 unsigned int cmp = *string_length-1;
 unsigned char previous;
 unsigned char counter;
 unsigned char*encoded = (unsigned char *)calloc(*string_length*2,1);
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
  if (counter>flag_runs)
  {
   //printf("%d,",counter);
   c=flag_runs;
   while (c>0)
   {
    encoded[j] = previous;
    j++;
    c--;
   }
   encoded[j] = counter-flag_runs-1;
   new_alphabet[k] = encoded[j];
   printf("-----%d----\n",counter);
   j++;
   k++;
  }
  else if (counter>1)
  {
    for (c=counter-1;c>0;c--)
    {
     encoded[j] = s[i-c];
     j++;
    }
  }
 }

 //printf("j = %d length = %d\n",j,*string_length);
 encoded[j]='\0';
 if (j!=*string_length)
  encoded = (unsigned char *)realloc(encoded,j+1);
 printf("usetrilo sa %d bytov\n",*string_length-j);
 *string_length = j;
 printf("hodnota string_length je %d\n",*string_length);
 new_alphabet[k]='\0';
 *new_alphabet_length = k;
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
 /*for (i=0;i<alphabet_length;i++)
 {
 	printf("%d=%d\n",heap[i]->freq,heap[i]->symbol);
 }
 printf("\n");*/
 for (i=alphabet_length/2;i>=0;i--)
 {
 	heapify(heap,i,alphabet_length);
 }
/*for (i=0;i<alphabet_length;i++)
 {
 	printf("%d=%d\n",heap[i]->freq,heap[i]->symbol);
 }
 printf("\n");*/

 return heap;
}

struct huffman_node *build_huffman_tree(unsigned char*alphabet,unsigned int*freq, unsigned char alphabet_length)
{
 struct huffman_node ** heap = build_min_heap(alphabet,freq,alphabet_length);
 struct huffman_node *left_node;
 struct huffman_node *right_node;
 struct huffman_node *root_node;
 struct huffman_node *current;
 unsigned int i = 0;
 unsigned int j = 0;
 unsigned int k;
 int size = alphabet_length;
 
 unsigned int huffman_tree_size = size+size-1;
 while (size>1)
 {
 //extract 2 min frequencies from heap
 left_node = get_root(heap,size);
 //printf("extrahovali sme %d %c\n",left_node->freq,left_node->symbol);
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
//printf("extrahovali sme %d %c\n",right_node->freq,right_node->symbol);
 //add extracted frequencies and its sum to huffman tree
 //add to heap sum of extracted roots of minheap
 root_node = create_new_node('a',left_node->freq+right_node->freq);
 root_node->left = left_node;
 root_node->right = right_node;
 heap[size] = root_node;
 size++;
 }
 return root_node;
}

//encode input string according to huffman tree
unsigned char* pack_huffman_to_bitvector(struct huffman_node *root_node, unsigned char *alphabet, unsigned char alphabet_length, unsigned char *s, unsigned int string_length, unsigned int*bitvector_length, unsigned char flag_new_alphabet)
{
 struct huffman_node *current = root_node;
 unsigned int i = 0;
 unsigned int j = 0;
 unsigned int k;
 unsigned int m;
 unsigned char code[alphabet_length];
 unsigned char result_codes[alphabet_length][alphabet_length];
 printf("idem vyhladat %d\n",current->symbol);
    fflush(stdout);
 unsigned char *bitvector = (unsigned char*)calloc(string_length,1);
 printf("%d %d\n",string_length,alphabet_length);
 fflush(stdout);
 unsigned int bits_counter = 0;
 unsigned int byte_index;
 unsigned char alphabet_index;

 //keep traversing huffman tree until are generated codes for all leaves 
 while (1)
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
  
  //if its leaf, insert code according to traversal from root to leaf
   if (current->left==NULL && current->right==NULL)
   {
    printf("idem vyhladat %d\n",current->symbol);
    fflush(stdout);
    i = get_index_in_alphabet(alphabet,current->symbol);
    printf("najdeny %d\n",i);
    fflush(stdout);
    for(k=0;k<j;k++)
    {
    result_codes[i][k] = code[k];
    }
    result_codes[i][k]='\0';
   }
   if (current==root_node)
    break;
 }

 //when all codes are generated for alphabet, encode each character of input string according to code
 j = 0;
 byte_index = 0;
 m = 0;
 for(i=0;i<string_length;i++)
 {
  if (flag_new_alphabet)
   alphabet_index = get_index_in_alphabet(alphabet,s[i]);
  else
   alphabet_index = s[i];
  bits_counter = j + strlen(result_codes[alphabet_index]);
  for (k=0;j<bits_counter;j++,k++,m++)
  {
    //if eight bit, go to next byte, otherwise shift left current byte for next value
   if (m==8)
   {
    byte_index++;
    m = 0;
   }
   else
    bitvector[byte_index] = bitvector[byte_index]<<1;
   bitvector[byte_index] = bitvector[byte_index] | (result_codes[alphabet_index][k]-'0');
   //printf("stav je %d, pripocitavam %d, vysledok je %d\n",bitvector[byte_index],result_codes[alphabet_index][k]-'0',bitvector[byte_index] | (result_codes[alphabet_index][k]-'0'));
   //printf("ukladam na poz %d bit %d, hodnotu %d\n",byte_index,m,result_codes[alphabet_index][k]);
  }
 }
 bitvector[byte_index] = bitvector[byte_index]<<(8-m);
 byte_index++;
 bitvector = (unsigned char *)realloc(bitvector,byte_index);
 printf("pocet bytov je %d\n",byte_index);
 *bitvector_length = j;
 return bitvector;
}

unsigned char* order_new_alphabet(unsigned char*new_alphabet,unsigned char*alphabet,unsigned char *result_length, unsigned int new_alphabet_length)
{
 unsigned int alphabet_length = strlen(alphabet);
 unsigned char *temp = (unsigned char*) malloc(new_alphabet_length);
 unsigned char c;
 unsigned char *result;
 int i;
 int j=0;
 unsigned int k;
 printf("alphabet_length je %d\n",alphabet_length);
 printf("new alphabet_length je %d\n",new_alphabet_length);
 for(i=0;i<new_alphabet_length;i++)
  printf("%d ",new_alphabet[i]);
printf("\n");
 for (i=0;i<new_alphabet_length;i++)
 {
  if (new_alphabet[i]>=alphabet_length)
  {
   if (j==0)
     temp[j++] = new_alphabet[i];
   else
   {
    k = 0;
    printf("new_alphabet[i] je %d, j je %d\n",new_alphabet[i],j);
    while (k<=j)
    {
     printf("k je %d, j je %d\n",new_alphabet[i],j);
     if (temp[k]==new_alphabet[i])
      break;
     k++;
    }
    if (temp[k]!=new_alphabet[i])
    {
     temp[j] = new_alphabet[i];
     j++;
    }
   }
  }
 }
 
 //order temp
 for (i = 0; i < j-1; i++) {
      for (k = i+1; k < j; k++) {
         if (temp[i] > temp[k]) {
            c = temp[i];
            temp[i] = temp[k];
            temp[k] = c;
         }
      }
   }
 
 result = (unsigned char*)malloc(alphabet_length+j);
 printf("velkost vyslednej abecedy je %d %d teda %ld\n",alphabet_length,j,strlen(result));
 for(i=0;i<alphabet_length;i++)
  result[i]=i;
 strncpy(&result[alphabet_length],temp,j);
 
 for(i=0;i<alphabet_length+j;i++)
  printf("%d ",result[i]);
 printf("\n");
  result[alphabet_length+j]='\0';
 *result_length = alphabet_length+j;
 free(temp);
 free(new_alphabet);
 return result;
}

void free_huffman_tree(struct huffman_node *node)
{
 if (node!=NULL)
 {
  free_huffman_tree(node->left);
  free_huffman_tree(node->right);
  free(node);
 }
}

struct wavelet_tree *build_huffman_shaped_WT(unsigned char *s, unsigned int *frequencies)
{
 //build huffman tree of alphabet
 struct huffman_node*root = build_huffman_tree(alphabet,frequencies,strlen(alphabet));

 //build wavelet tree recursively from root
 struct wavelet_tree *wt_root = build_WT_node(root,s);
 
 return wt_root;
}

//build wavelet tree node and its children
struct wavelet_tree *build_WT_node(struct huffman_node *root, unsigned char *s)
{
 if (root->left == NULL && root->right==NULL)
    return NULL;

 unsigned char * left_alphabet = get_alphabet(root->left);
 unsigned char * right_alphabet = get_alphabet(root->right);
 unsigned char max_bits = sizeof(unsigned long long int)*8;

 unsigned long long int *bitvector = (unsigned long long int*)calloc(genome_length/max_bits+1,sizeof(unsigned long long int*));
 //printf("size je %d a pocet slov je %d, bajtov je %d\n",string_length,string_length/max_bits+1,(string_length/max_bits+1)*8);
 
 unsigned int word_index = 0;
 unsigned char bit_index = 0;
 unsigned int i;

 printf("coding as 1: %s, ",right_alphabet);
 printf("coding as 0: %s, ",left_alphabet);
 for (i=0;i<genome_length;i++)
 {
  if (in_alphabet(s[i],right_alphabet))
  {
    //zakoduj ako 1
    //printf("1");
    bitvector[word_index] = bitvector[word_index]<<1;
    bitvector[word_index] = bitvector[word_index] + 1;
    bit_index++;
  }
  else if (in_alphabet(s[i],left_alphabet))
  {
    //zakoduj ako 0
    //printf("0");
    bitvector[word_index] = bitvector[word_index]<<1;
    bit_index++;
  }
  if (bit_index==max_bits){
    bit_index = 0;
    //printf("\n%llu\n",bitvector[word_index]);
    word_index++;

  }
 }
 //printf("pocet slov alokovanych je %d, i je %d\n",word_index,i);
 bitvector[word_index] = bitvector[word_index]<<(max_bits-bit_index);
 //printf("idem realokovat na %ld\n",(word_index+1)*sizeof(unsigned long long int));
 //fflush(stdout);
 bitvector = realloc(bitvector, (word_index+1)*sizeof(unsigned long long int));
 printf(" - used %ld bytes\n",(word_index+1)*sizeof(unsigned long long int));
 struct wavelet_tree* wt_node = (struct wavelet_tree*)malloc(sizeof(struct wavelet_tree));
 wt_node->right_alphabet = right_alphabet;
 wt_node->left_alphabet = left_alphabet;
 wt_node->bitvector = bitvector;
 if (root->left!=NULL)
  wt_node->left = build_WT_node(root->left,s);
 else
  wt_node->left = NULL;
 if (root->right!=NULL)
  wt_node->right = build_WT_node(root->right,s);
 else
  wt_node->right = NULL;

wt_node->bitcount_table = build_bitcount_table(wt_node->bitvector,word_index+1);
return wt_node;
}

//build table of counters of bits in bitvectors of wavelet tree
unsigned int* build_bitcount_table(unsigned long long int *bitvector, unsigned int bitvector_length)
{
 unsigned int i;
 unsigned int j=0;
 unsigned int count = 0;
 unsigned int number = 1+(bitvector_length-1)/sample_OCC_size;
 unsigned int* bitcount_table = (unsigned int*) malloc (sizeof(unsigned int)*number);
 printf("building bitcount table on %d words(64bits) - used %d bytes\n",bitvector_length,number*sizeof(unsigned int));
 for (i = 0; i<bitvector_length;i++)
 {
  count = count + __builtin_popcountll(bitvector[i]);
  if ((i%sample_OCC_size)==(sample_OCC_size-1))
  {
   fflush(stdout);
   bitcount_table[j++] = count;
  }
 }
 return bitcount_table;
}

//function which returns alphabet which is coded in selected node of wavelet tree
unsigned char *get_alphabet(struct huffman_node*node)
{
 unsigned char *alphabet = NULL;
 unsigned char *left_alphabet = NULL;
 unsigned char *right_alphabet = NULL;
 if(node->left==NULL && node->right==NULL)
 {
  fflush(stdout);
  alphabet = (unsigned char *)malloc(8);
  if (alphabet==NULL){
    printf("error pri alokovani\n");
    fflush(stdout);
  }
  alphabet[0] = node->symbol;
  alphabet[1] = '\0';
 }
 else
 {
    left_alphabet = get_alphabet(node->left);
    right_alphabet = get_alphabet(node->right);
    const size_t left_alphabet_size = strlen(left_alphabet);
    const size_t right_alphabet_size = strlen(right_alphabet);

    alphabet = (unsigned char *)malloc(left_alphabet_size+right_alphabet_size+2);

    memcpy(alphabet, left_alphabet, left_alphabet_size);
  memcpy(alphabet+left_alphabet_size, right_alphabet, right_alphabet_size+1);
  free(left_alphabet);
  free(right_alphabet);
 }
 return alphabet;
}

//help function which checks if input character is in input alphabet
unsigned char in_alphabet(unsigned char c, unsigned char *alphabet)
{
 unsigned char i;
 for (i=0;i<strlen(alphabet);i++){
  if (alphabet[i]==c)
    return 1;
 }
 return 0;
}

//function which returns character on input position
unsigned char wt_access(unsigned int position, struct wavelet_tree *root)
{
 unsigned int word_index = position/max_bits;
 unsigned int bit_index = position-word_index*max_bits;
//if bit of input position is 1
 if (get_bit(root->bitvector[word_index],bit_index))
 {
  //printf("pristupujem R pos %d, %d %d %llu, %d \n",position,word_index,bit_index,current->bitvector[word_index],get_bit(current->bitvector[word_index],bit_index));
  //check if current node is leaf, therefore coded alphabet is only one character
   if (strlen(root->right_alphabet)==1)
  {
    return root->right_alphabet[0];
  }
  //if not in leaf, check right child, because bit = 1, and new position is count of 1
  return wt_access(count_set_bits(root->bitvector,position,root->bitcount_table)-1,root->right);
 }
 //vpravo CG C-10 G-11
 //vlaov A   A-00 T-01

 //if bit of input position is 0
 else
 {
  //printf("pristupujem L pos%d, %d %d %llu, %d \n",position,word_index,bit_index,current->bitvector[word_index],get_bit(current->bitvector[word_index],bit_index));
  
  //check if current node is leaf, therefore coded alphabet is only one character
  if (strlen(root->left_alphabet)==1)
  {
   return root->left_alphabet[0];
  }
  //if not in leaf, check left child, because bit = 1, and new position is count of 0
  return wt_access(count_unset_bits(root->bitvector,position,root->bitcount_table)-1,root->left);
 }
}

//function which counts 1 in bitvector up to input position while using auxiliary srtuctrre of bitcount table
unsigned int count_set_bits(unsigned long long int*bitvector, unsigned int position, unsigned int* bitcount_table)
{
 unsigned long long int remainder;
 unsigned int count = 0;
 unsigned int i = 0;

 //get index of nearest stored counter of ones in bitvector
 int index = position/(sample_OCC_size*max_bits)-1;
 //printf("index je %d\n",index);

 if (index>=0)
 {
  count = bitcount_table[index];
  i = (index+1)*sample_OCC_size;
  position -= i*max_bits;
 }

 //count 1 in remaining part of bitvector
 while (position>=max_bits)
 {
  count += __builtin_popcountll(bitvector[i++]);
  position = position - max_bits;
 }
 //printf("bitvector je %llu\n",bitvector[i]);
 //fflush(stdout);
 remainder = bitvector[i] >> (max_bits - position - 1);
 //printf("count je %d\n",count);
 //fflush(stdout);
 count+= __builtin_popcountll(remainder);
 //printf("count je %d\n",count);
return count;
}

unsigned int count_unset_bits(unsigned long long int*bitvector, unsigned int position, unsigned int* bitcount_table)
{
 unsigned long long int remainder;
 unsigned int count = 0;
 unsigned int i = 0;

 //printf("UNS position je %d,max_bits je %d, sample_occ_size je %d\n",position,max_bits,sample_occ_size);
 int index = position/(sample_OCC_size*max_bits)-1;
 //printf("index je %d, i je %d occ je %d\n",index,i,bitcount_table[index]);
 if (index>=0)
 {
  i = (index+1)*sample_OCC_size;
  count = i*max_bits - bitcount_table[index];
  position = position - i*max_bits;
 }

 //printf("index je %d, i je %d tvysnok je %d\n",index,i,position);
 
 while (position>=max_bits)
  {
   count = count + max_bits - __builtin_popcountll(bitvector[i++]);
   position = position - max_bits;
  }
  //printf("bitvector je %llu\n",bitvector[i]);
  remainder = bitvector[i] >> (max_bits-position - 1);
  //printf("count je %d\n",count);
  count += position + 1 - __builtin_popcountll(remainder);
  //printf("count je %d\n",count);
 return count;
}

unsigned long long int get_bit(unsigned long long int var, unsigned int position)
{
//return ((var) & (1LL<<(64-position)));
  return (var >> (63-position)) & 1;
}

unsigned int wt_rank(unsigned char c, unsigned int position, struct wavelet_tree *root)
{
 unsigned int count = 0;
 if (root == NULL)
 {
  return position;
 }

 if (position)
    --position;
else return 0;
 //printf("WTF %c, position je %d\n",c,position);
//fflush(stdout);
 //ak sa znak nachadza v lavej abecede spocitaju sa nuly inac 1
 if (in_alphabet(c,root->left_alphabet))
 {
  //printf("pocitam 0 %c do pozicie %d\n",c,position);
  count = count_unset_bits(root->bitvector,position,root->bitcount_table);
  return wt_rank(c,count,root->left);
 }
 else
 {
  //printf("pocitam 1 %c do pozicie %d\n",c,position);
  count = count_set_bits(root->bitvector,position,root->bitcount_table);
  return wt_rank(c,count,root->right);
 }
}