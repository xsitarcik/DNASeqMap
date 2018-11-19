#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "wavelet_tree.h"

struct wavelet_tree *build_huffman_shaped_WT(unsigned char *s, unsigned int *frequencies)
{
 //build huffman tree of alphabet
 struct huffman_node*root = build_huffman_tree((unsigned char *)alphabet,frequencies,strlen(alphabet));

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
  if (in_alphabet(s[i],(char *)right_alphabet))
  {
    //zakoduj ako 1
    //printf("1");
    bitvector[word_index] = bitvector[word_index]<<1;
    bitvector[word_index] = bitvector[word_index] + 1;
    bit_index++;
  }
  else if (in_alphabet(s[i],(char *)left_alphabet))
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

unsigned long long int get_bit(unsigned long long int var, unsigned int position)
{
//return ((var) & (1LL<<(64-position)));
  return (var >> (63-position)) & 1;
}

unsigned int get_bit_4byte(unsigned int var, unsigned int position)
{
//return ((var) & (1LL<<(64-position)));
  return (var >> (31-position)) & 1;
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
 if (in_alphabet(c,(char *)root->left_alphabet))
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
   if (strlen((char *)root->right_alphabet)==1)
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
  if (strlen((char *)root->left_alphabet)==1)
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
   count += __builtin_popcountll(~bitvector[i++]);
   position = position - max_bits;
  }
  //printf("bitvector je %llu\n",bitvector[i]);
  remainder = bitvector[i] >> (max_bits-position - 1);
  //printf("count je %d\n",count);
  count += position + 1 - __builtin_popcountll(remainder);
  //printf("count je %d\n",count);
 return count;
}

//build table of counters of bits in bitvectors of wavelet tree
unsigned int* build_bitcount_table(unsigned long long int *bitvector, unsigned int bitvector_length)
{
 unsigned int i;
 unsigned int j=0;
 unsigned int count = 0;
 unsigned int number = 1+(bitvector_length-1)/sample_OCC_size;
 unsigned int* bitcount_table = (unsigned int*) malloc (sizeof(unsigned int)*number);
 printf("building bitcount table on %d words(64bits) - used %lu bytes\n",bitvector_length,number*sizeof(unsigned int));
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

unsigned char get_set_bits(unsigned char bits)
{
 int ret = 1;
 for (unsigned int i=0;i<bits;i++)
  ret = ret * 2;
 return (ret - 1);
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


//help function which checks if input character is in input alphabet
unsigned char in_alphabet(unsigned char c, char *alphabet)
{
 unsigned char i;
 for (i=0;i<strlen(alphabet);i++){
  if (alphabet[i]==c)
    return 1;
 }
 return 0;
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
    const size_t left_alphabet_size = strlen((char *)left_alphabet);
    const size_t right_alphabet_size = strlen((char *)right_alphabet);

   alphabet = (unsigned char *)malloc(left_alphabet_size+right_alphabet_size+2);

  memcpy(alphabet, left_alphabet, left_alphabet_size);
  memcpy(alphabet+left_alphabet_size, right_alphabet, right_alphabet_size+1);
  free(left_alphabet);
  free(right_alphabet);
 }
 return alphabet;
}
