struct symbol_table
{
 char symbol;
 struct symbol_table *next;
 struct symbol_table *previous;
};

struct symbol_table *build_symbol_table(char *alphabet);
struct symbol_table *push_to_front(struct symbol_table *front, struct symbol_table *current);
unsigned char get_min_bits_per_char(char *alphabet);
void move_to_front_encode (char *alphabet, char *s);
char *move_to_front_decode (char *alphabet, unsigned int bitvector_length,unsigned char*bitvector);
unsigned char get_set_bits(unsigned char bits);
void print_bit_vector(unsigned char *bitvector, unsigned int bitvector_length);
struct symbol_table *decode(struct symbol_table *front, unsigned char code);
unsigned char *arithmetic_encode (unsigned int *bitvector_length, char *s, char * alphabet, unsigned int string_length);
unsigned char* zero_runs_encode(unsigned char *s, unsigned int *string_length);
