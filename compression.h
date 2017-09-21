struct symbol_table
{
 char symbol;
 struct symbol_table *next;
 struct symbol_table *previous;
};

struct symbol_table *build_symbol_table(char *alphabet);
unsigned char *move_to_front_encode (struct symbol_table *front, char *s, unsigned int *bitvector_length);
unsigned char get_set_bits(unsigned char bits);
void print_bit_vector(unsigned char *bitvector, unsigned int bitvector_length);
