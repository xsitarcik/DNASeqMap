extern char *alphabet;
extern unsigned char alphabet_size;
extern unsigned char max_error;
extern unsigned int sample_OCC_size; //in reality it's *64
extern unsigned int sample_SA_size;
extern unsigned int genome_length;
extern unsigned char max_bits;
extern struct FMIndex_WT *FM_index_WT;
extern unsigned char*genome;
extern struct wavelet_tree *WT_root;
extern unsigned int*count_table;
extern unsigned int*entries;
extern unsigned int THRESHOLD;

struct FMIndex
{
 unsigned int bitvector_length;
 unsigned char *alphabet;
 unsigned int end;
 unsigned int *sampleSA;
 unsigned int *count_table;
 unsigned int **occurence_table;
 unsigned char *bwt;
 unsigned char alphabetically_encoded;
};

struct FMIndex_WT
{
unsigned int end;
unsigned int *sampleSA;
};

struct compressedFMIndex
{
 unsigned int length;
 unsigned int block_size;
 unsigned char *alphabet;
 unsigned char flag_mtf;
 unsigned char flag_runs;
 unsigned char flag_huffman;
 unsigned int end;
 unsigned int *sampleSA;
 unsigned int *count_table;
 struct compressed_block *array_of_blocks;
};


long long int approximate_search_in_FM_index_entry(unsigned char *pattern, unsigned int*result);
void threshold_search_pattern_in_FM_index_entry(char *pattern, unsigned int *result_length, unsigned int *result, unsigned int pattern_length);
unsigned int extend_seed_in_FM_index_entry(unsigned char*pattern, unsigned int last, unsigned int*result);
unsigned char search_pattern_in_FM_index_entry(char *pattern, unsigned char current_pattern_length, unsigned int*result);
unsigned int get_SA_value_entry(unsigned int bwt_position);
unsigned char wt_access_entry(unsigned int entry_index, unsigned short int in_entry_index);
unsigned int wt_rank_entry(unsigned char c, unsigned int position, unsigned int entry_index, unsigned int in_entry_index);
void rebuild_FM_index_into_entries(unsigned int*suffix_array, unsigned char*bwt);
struct FMIndex_WT*build_FM_index_WT(unsigned int *suffix_array, unsigned char *bwt);
struct FMIndex*build_FM_index(unsigned int *suffix_array, unsigned char *bwt);
struct compressedFMIndex*build_compressed_FM_index(unsigned int *suffix_array, unsigned char *bwt,
  unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman, unsigned int block_size);
unsigned int find_end(unsigned int *suffix_array);
unsigned int *create_sample_SA(unsigned int *suffix_array);
unsigned int *create_count_table(unsigned char *s);
unsigned int get_SA_value(unsigned int bwt_position, unsigned char c, struct FMIndex *fm_index);
unsigned int get_SA_value_compressed(unsigned int bwt_position, unsigned int length, unsigned char *alphabet, unsigned int*sampleSA, 
  unsigned int *count_table,unsigned int block_size, struct compressed_block*block, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman);
unsigned int **create_occurence_table(unsigned char *s, unsigned char *alphabet);
unsigned int count_occ(unsigned char *s, unsigned int **occurence_table, unsigned int position, unsigned char c, unsigned char character);
unsigned int last_to_first(unsigned char c, unsigned int bwt_position, unsigned char *alphabet,unsigned int *count_table,unsigned char*bwt,unsigned int**occurence_table);
unsigned char get_alphabet_index(char *alphabet, unsigned char c);
unsigned char *reverseBWT(unsigned int end, unsigned char*alphabet,unsigned int*count_table,unsigned char*bwt, unsigned int**occurence_table,unsigned char alphabetically_encoded);
void print_occurence_table(struct FMIndex *fm_index);
void print_sample_SA(struct FMIndex *fm_index);
void print_count_table(struct FMIndex *fm_index);
void print_info_fm_index(struct FMIndex *fm_index);
unsigned int*search_pattern(struct FMIndex *fm_index, char *pattern);
unsigned int*approximate_search(struct FMIndex *fm_index, unsigned char *pattern);
unsigned int*approximate_search_in_compressed_FM_index(struct compressedFMIndex *compressed_fm_index, unsigned char *pattern, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman);
unsigned int score(char a, char b);
unsigned char align(char *p1, char*p2);
unsigned int calculate(int i, int j, int k);
unsigned int get_max_array(unsigned char* array, unsigned int length);
unsigned int last_to_first_encoded(unsigned char c, unsigned int bwt_position, unsigned char *alphabet,unsigned int *count_table,unsigned char*bwt,unsigned int**occurence_table);
unsigned int last_to_first_in_compressed_FMIndex(unsigned char c, unsigned int bwt_position, unsigned char*alphabet, unsigned int* count_table,struct compressed_block *block, unsigned int block_size, 
  unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman);
unsigned char *decompress_FMIndex(struct compressedFMIndex *compressed_FM_index);
unsigned int count_occ_in_block(struct compressed_block *block, unsigned int position,unsigned char c,unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman, unsigned char*alphabet);
unsigned int count_occ_in_compressed_FMIndex(struct compressed_block *block, unsigned int block_size, unsigned int position, unsigned char c, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman, unsigned char*alphabet);
unsigned int count_occ_in_decompressed_FMIndex(struct compressed_block *block, unsigned char*bwt,unsigned int block_size, unsigned int position, unsigned char c, unsigned char*alphabet);
unsigned int*search_pattern_in_compressed_FM_index(struct compressedFMIndex *compressed_fm_index, char *pattern,unsigned char flag_mtf,unsigned char flag_runs, unsigned char flag_huffman);
unsigned char search_pattern_in_FM_index_WT(char *pattern, unsigned char current_pattern_length, unsigned int*result);
long long int approximate_search_in_FM_index_WT(unsigned char *pattern, unsigned int*result);
unsigned int get_SA_value_WT(unsigned int bwt_position);
unsigned int last_to_first_WT(unsigned char c, unsigned int bwt_position, struct wavelet_tree *wtree_root,unsigned int *count_table);
void threshold_search_pattern_in_FM_index_WT(char *pattern, unsigned int *result_length, unsigned int *result,unsigned int pattern_length);
unsigned int extend_seed_in_FM_index_WT(unsigned char*pattern, unsigned int last, unsigned int*result);