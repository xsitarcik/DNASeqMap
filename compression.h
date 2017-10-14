struct symbol_table
{
 char symbol;
 struct symbol_table *next;
 struct symbol_table *previous;
};

unsigned char*  compress (unsigned char flag_mtf, unsigned char flag_zero_runs,unsigned char huffman, unsigned char*alphabet, unsigned char*bwt, unsigned int*length, unsigned char *alphabetically_encoded);
struct symbol_table *build_symbol_table(char *alphabet);
struct symbol_table *push_to_front(struct symbol_table *front, struct symbol_table *current);
unsigned char get_min_bits_per_char(char *alphabet);
void move_to_front_encode (char *alphabet, char *s);
char *move_to_front_decode (char *alphabet, unsigned int bitvector_length,unsigned char*bitvector);
unsigned char get_set_bits(unsigned char bits);
void print_bit_vector(unsigned char *bitvector, unsigned int bitvector_length);
struct symbol_table *decode(struct symbol_table *front, unsigned char code);
unsigned char *arithmetic_encode (char *s, char * alphabet, unsigned int *string_length);
unsigned char* zero_runs_encode(unsigned char *s, unsigned int *string_length);
unsigned char get_index_in_alphabet(char *alphabet, char c);
void alphabet_encode (char *s, char * alphabet);
unsigned char *arithmetic_decode (unsigned char *bitvector, char *alphabet, unsigned int *bitvector_length);
unsigned char decode_bits(unsigned char bits_per_char,unsigned int bitposition, unsigned char *bitvector);
