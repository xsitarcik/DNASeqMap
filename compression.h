struct symbol_table
{
 unsigned int index;
 char symbol;
 struct symbol_table *next;
 struct symbol_table *previous;
};

struct symbol_table *build_symbol_table(char *alphabet);
unsigned int move_to_front_encode (struct symbol_table *front, char *s);
