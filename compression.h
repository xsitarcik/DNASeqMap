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
 unsigned int*runs;
};

struct huffman_node
{
 unsigned int freq;
 unsigned char symbol;
 struct huffman_node*left;
 struct huffman_node*right;
};


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
unsigned char *arithmetic_encode (char *s, unsigned int string_length, unsigned char bits_per_char, unsigned int *bitvector_length);
unsigned char* run_length_encode(unsigned char *s, unsigned int *string_length);
unsigned char get_index_in_alphabet(char *alphabet, char c);
void alphabet_encode (char *s, char * alphabet,unsigned int block_size);
unsigned char *arithmetic_decode (unsigned char *bitvector, char *alphabet, unsigned int string_length);
unsigned char decode_bits(unsigned char bits_per_char,unsigned int bitposition, unsigned char *bitvector);
unsigned int *calc_occurences(char *s, int string_length, char *alphabet);
unsigned char *decompress_block(unsigned int bitvector_length, unsigned char*bitvector, unsigned int*runs, unsigned char flag_mtf, 
	unsigned char flag_runs, unsigned char flag_huffman, unsigned int block_size, unsigned char*alphabet);
unsigned char* reverseBWT_compressed(unsigned char*bwt, unsigned int length, unsigned int end, unsigned int* count_table, unsigned char*alphabet, 
  unsigned int block_size, struct compressed_block *block, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman);
struct huffman_node **build_min_heap(unsigned char*alphabet, unsigned int*frequencies,unsigned int alphabet_length);
void heapify(struct huffman_node **heap, unsigned int index, unsigned int size);
struct huffman_node *build_huffman_tree(unsigned char*alphabet,unsigned int*freq);