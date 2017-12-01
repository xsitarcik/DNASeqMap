unsigned int *init_suffix_array(unsigned int *suffix_array,unsigned char *s, unsigned int genome_length);
unsigned int max(unsigned int number1,unsigned int number2);
unsigned int min(unsigned int number1,unsigned int number2);
unsigned int compare_rotations(unsigned char *s, unsigned int start1, unsigned int start2, unsigned int genome_length);
unsigned int *insertion_sort_array(unsigned int *suffix_array, unsigned char *s, unsigned int genome_length);
unsigned char *create_bwt(unsigned int *suffix_array, unsigned char *s, unsigned int genome_length);
unsigned char*load_genome_from_file(unsigned char*file,unsigned int*genome_length);
void reverse_string(unsigned char *str);