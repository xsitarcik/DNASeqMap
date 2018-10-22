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

struct wavelet_tree *build_WT_node(struct huffman_node *root, unsigned char *s);
unsigned long long int get_bit(unsigned long long int var, unsigned int position);
unsigned int get_bit_4byte(unsigned int var, unsigned int position);
unsigned int wt_rank(unsigned char c, unsigned int position, struct wavelet_tree *root);
unsigned char wt_access(unsigned int position, struct wavelet_tree *root);
unsigned int count_set_bits(unsigned long long int*bitvector, unsigned int position, unsigned int* bitcount_table);
unsigned int count_unset_bits(unsigned long long int*bitvector, unsigned int position, unsigned int* bitcount_table);
unsigned int* build_bitcount_table(unsigned long long int *bitvector, unsigned int bitvector_length);
unsigned char get_set_bits(unsigned char bits);
unsigned char *bit_pack (char *s, unsigned int string_length, unsigned char bits_per_char, unsigned int *bitvector_length);
unsigned char *get_alphabet(struct huffman_node*node);
unsigned char in_alphabet(unsigned char c, unsigned char *alphabet);


