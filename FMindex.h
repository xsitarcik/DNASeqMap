extern unsigned char *alphabet;
extern unsigned char alphabet_size;
extern unsigned char max_error;
extern unsigned int sample_OCC_size; //in reality it's *64
extern unsigned int sample_SA_size;
extern unsigned int genome_length;
extern unsigned char max_bits;

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
unsigned int *count_table;
struct wavelet_tree *WT_root;
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
unsigned int get_alphabet_index(unsigned char *alphabet, unsigned char c);
unsigned char *reverseBWT(unsigned int end, unsigned char*alphabet,unsigned int*count_table,unsigned char*bwt, unsigned int**occurence_table,unsigned char alphabetically_encoded);
void print_occurence_table(struct FMIndex *fm_index);
void print_sample_SA(struct FMIndex *fm_index);
void print_count_table(struct FMIndex *fm_index);
void print_info_fm_index(struct FMIndex *fm_index);
unsigned int*search_pattern(struct FMIndex *fm_index, char *pattern);
unsigned int*approximate_search(struct FMIndex *fm_index, unsigned char *pattern);
unsigned int*approximate_search_in_compressed_FM_index(struct compressedFMIndex *compressed_fm_index, unsigned char *pattern, unsigned char flag_mtf, unsigned char flag_runs, unsigned char flag_huffman);
unsigned int score(char a, char b);
void align(char *p1, char*p2, int error);
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
unsigned int*search_pattern_in_FM_index_WT(struct FMIndex_WT *FM_index_WT, char *pattern);
long long int approximate_search_in_FM_index_WT(struct FMIndex_WT *FM_index_WT, unsigned char *pattern, unsigned int *result);
unsigned int get_SA_value_WT(unsigned int bwt_position, struct FMIndex_WT *fm_index_wt);
unsigned int last_to_first_WT(unsigned char c, unsigned int bwt_position, struct wavelet_tree *wtree_root,unsigned int *count_table);
unsigned int*threshold_search_pattern_in_FM_index_WT(struct FMIndex_WT *FM_index_WT, char *pattern, unsigned int *result_length, unsigned int *result);