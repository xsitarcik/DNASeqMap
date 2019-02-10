extern char *alphabet;
extern unsigned char alphabet_size;
extern unsigned int sample_OCC_size; //in reality it's *64
extern unsigned int sample_SA_size;
extern unsigned int genome_length;
extern unsigned char max_bits;

struct huffman_node
{
 unsigned char visited;
 unsigned short int freq;
 unsigned char symbol;
 struct huffman_node*left;
 struct huffman_node*right;
};

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))

struct wavelet_tree
{
 struct wavelet_tree *right;
 struct wavelet_tree *left;
 unsigned char *left_alphabet;
 unsigned char *right_alphabet;
 unsigned long long int *bitvector;
 unsigned int *bitcount_table;
};

struct wavelet_tree *build_huffman_shaped_WT(unsigned char *s, unsigned int *frequencies);
struct wavelet_tree *build_WT_node(struct huffman_node *root, unsigned char *s);
unsigned long long int get_bit(unsigned long long int var, unsigned int position);
unsigned int get_bit_4byte(unsigned int var, unsigned int position);
unsigned int get_bit_4byte_reversed(unsigned int var, unsigned int position);
unsigned int wt_rank(unsigned char c, unsigned int position, struct wavelet_tree *root);
unsigned char wt_access(unsigned int position, struct wavelet_tree *root);
unsigned int count_set_bits(unsigned long long int*bitvector, unsigned int position, unsigned int* bitcount_table);
unsigned int count_unset_bits(unsigned long long int*bitvector, unsigned int position, unsigned int* bitcount_table);
unsigned int* build_bitcount_table(unsigned long long int *bitvector, unsigned int bitvector_length);
unsigned char get_set_bits(unsigned char bits);
unsigned char *bit_pack (char *s, unsigned int string_length, unsigned char bits_per_char, unsigned int *bitvector_length);
unsigned char in_alphabet(unsigned char c, char *alphabet);
struct huffman_node *build_huffman_tree(unsigned char*alphabet, unsigned int*freq, unsigned char alphabet_length);
struct huffman_node **build_min_heap(unsigned char*alphabet, unsigned int*frequencies,unsigned int alphabet_length);
void heapify(struct huffman_node **heap, unsigned int index, unsigned int size);
struct huffman_node *build_huffman_tree(unsigned char*alphabet, unsigned int*freq, unsigned char alphabet_length);
struct huffman_node *create_new_node(unsigned char c, unsigned int freq);
struct huffman_node *get_root(struct huffman_node **heap, unsigned int size);
unsigned char *get_alphabet(struct huffman_node*node);