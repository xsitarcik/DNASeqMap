struct FMIndex
{
 int *sampleSA;
 int *count_table;
 int **occurence_table;
 char *bwt;
 int sample_SA_size;
 int sample_OCC_size;
};

int *create_sample_SA(int *suffix_array,int sample_size, int array_size);
int *create_count_table(char *s, int string_length, char* alphabet);
int get_SA_value(int bwt_position, char c, int character, struct FMIndex *fm_index);
int **create_occurence_table(char *s, int string_length, char *alphabet, int sample_size);
int count_occ(char *s, int **occurence_table, int position, char c, int character, int sample_size);
int last_to_first(char c, int character, int bwt_position, struct FMIndex *fm_index);
int get_alphabet_index(char *alphabet, char c);
char *reverseBWT(char *bwt,int end,int length, char *alphabet, struct FMIndex *fm_index);
void print_occurence_table(int **occurence_table, int alphabet_size, int sample_OCC_size, int length);
