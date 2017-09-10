struct FMIndex
{
 int length;
 int sample_SA_size;
 int sample_OCC_size;
 char *alphabet;
 int *sampleSA;
 int *count_table;
 int **occurence_table;
 char *bwt;
};

struct FMIndex*build_FM_index(int *suffix_array, int sample_SA_size, int sample_OCC_size, int genome_length, char *bwt, char *alphabet);
int *create_sample_SA(int *suffix_array,int sample_size, int array_size);
int *create_count_table(char *s, int string_length, char* alphabet);
int get_SA_value(int bwt_position, char c, struct FMIndex *fm_index);
int **create_occurence_table(char *s, int string_length, char *alphabet, int sample_size);
int count_occ(char *s, int **occurence_table, int position, char c, int character, int sample_size);
int last_to_first(char c, int bwt_position, struct FMIndex *fm_index);
int get_alphabet_index(char *alphabet, char c);
char *reverseBWT(int end, struct FMIndex *fm_index);
void print_occurence_table(int **occurence_table, int alphabet_size, int sample_OCC_size, int length);
void print_sample_SA(struct FMIndex *fm_index);
