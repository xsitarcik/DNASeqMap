extern char *alphabet;
extern unsigned char alphabet_size;
extern unsigned char max_error;
extern unsigned int sample_OCC_size; //in reality it's *64
extern unsigned int sample_SA_size;
extern unsigned int genome_length;
extern unsigned char max_bits;

struct symbol_table
{
 char symbol;
 struct symbol_table *next;
 struct symbol_table *previous;
};

struct compressed_block
{
 unsigned int block_size;
 unsigned int bitvector_length;
 unsigned int*occurences;
 unsigned char*bitvector;
 struct huffman_node*huffman_tree;
};

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
unsigned char *get_alphabet(struct huffman_node*node);
unsigned char in_alphabet(unsigned char c, unsigned char *alphabet);
unsigned int wt_rank(unsigned char c, unsigned int position, struct wavelet_tree *root);
unsigned char wt_access(unsigned int position, struct wavelet_tree *root);
unsigned int count_set_bits(unsigned long long int*bitvector, unsigned int position, unsigned int* bitcount_table);
unsigned int count_unset_bits(unsigned long long int*bitvector, unsigned int position, unsigned int* bitcount_table);
unsigned int* build_bitcount_table(unsigned long long int *bitvector, unsigned int bitvector_length);


struct compressed_block *compress_FMIndex(unsigned int block_size, unsigned char flag_mtf, unsigned char flag_runs,unsigned char flag_huffman, 
	unsigned char*alphabet, unsigned char*bwt, unsigned int*length);
struct symbol_table *build_symbol_table(char *alphabet);
struct symbol_table *push_to_front(struct symbol_table *front, struct symbol_table *current);
unsigned char get_min_bits_per_char(char *alphabet);
void *move_to_front_encode (char *alphabet, char *s, unsigned int block_size);
void *move_to_front_decode (char *alphabet, unsigned int string_length,unsigned char*s);
unsigned char get_set_bits(unsigned char bits);
void print_bit_vector(unsigned char *bitvector, unsigned int bitvector_length);
struct symbol_table *decode(struct symbol_table *front, unsigned char code);
unsigned char *bit_pack (char *s, unsigned int string_length, unsigned char bits_per_char, unsigned int *bitvector_length);
unsigned char* run_length_encode(unsigned char *s, unsigned int *string_length, unsigned char* new_alphabet, unsigned int *new_alphabet_length, unsigned char flag_runs);
unsigned char get_index_in_alphabet(char *alphabet, char c);
void alphabet_encode (char *s, char * alphabet,unsigned int block_size);
unsigned char *arithmetic_decode_with_RLE (unsigned char *bitvector, char *alphabet, unsigned int string_length);
unsigned char *arithmetic_decode_without_RLE (unsigned char *bitvector, char *alphabet, unsigned int string_length);
unsigned char decode_bits(unsigned char bits_per_char,unsigned int bitposition, unsigned char *bitvector);
unsigned int *calc_occurences(char *s, int string_length, char *alphabet);
unsigned char *decompress_block(unsigned int bitvector_length, unsigned char*bitvector, unsigned char flag_mtf, 
	unsigned char flag_runs, unsigned char flag_huffman, unsigned int block_size, unsigned char*alphabet, struct huffman_node*huffman_tree);
unsigned char* reverseBWT_compressed(unsigned char*bwt, unsigned int length, unsigned int end, unsigned int* count_table, unsigned char*alphabet, 
  unsigned int block_size, struct compressed_block *block, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman);
struct huffman_node **build_min_heap(unsigned char*alphabet, unsigned int*frequencies,unsigned int alphabet_length);
void heapify(struct huffman_node **heap, unsigned int index, unsigned int size);
struct huffman_node *build_huffman_tree(unsigned char*alphabet, unsigned int*freq, unsigned char alphabet_length);
unsigned char* pack_huffman_to_bitvector(struct huffman_node *root_node, unsigned char *alphabet,unsigned char alphabet_length, unsigned char *s, unsigned int string_length, unsigned int*bitvector_length,unsigned char flag_new_alphabet);
unsigned char* order_new_alphabet(unsigned char*new_alphabet,unsigned char*alphabet,unsigned char *result_length, unsigned int new_alphabet_length);
unsigned char* huffman_decode_without_RLE(unsigned char*bitvector, struct huffman_node*huffman_tree,unsigned int string_length);
unsigned char* huffman_decode_with_RLE(unsigned char*bitvector, struct huffman_node*huffman_tree,unsigned int string_length,unsigned char flag_runs);
void free_huffman_tree(struct huffman_node *node);
unsigned char * naive_compress(unsigned char *genome);