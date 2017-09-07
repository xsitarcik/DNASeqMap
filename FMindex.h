struct FMIndex
{
 int *sampleSA;
 int *count_table;
 int **occurence_table;
 char *bwt;
};

int *create_sample_SA(int *suffix_array,int sample_size, int array_size);
int *create_count_table(char *s, int string_length, char* alphabet);
int get_SA_value(int *sample_SA, int sample_size, int bwt_position);
int **create_occurence_table(char *s, int string_length, char *alphabet, int sample_size);
int count_occ(char *s, int **occurence_table, int position, char c, int character, int sample_size);
int last_to_first(int bwt_position);
