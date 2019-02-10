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
extern unsigned int k_mers_permutation;
extern unsigned int total_kmers;
extern unsigned int* kmers_hash;
extern unsigned int*root_rankvector;
extern unsigned int*left_rankvector;
extern unsigned int*right_rankvector;
extern unsigned int*sample_SA;
/*extern clock_t t_get_sa_value; 
extern clock_t t_rank; 
extern clock_t t_access; 
extern clock_t t_approx_search;
extern clock_t t_align;  */

struct FMIndex_WT
{
unsigned int end;
unsigned int *sampleSA;
};

unsigned int get_SA_value_exSA(unsigned int bwt_position);
unsigned char search_pattern_in_FM_index_exSA(char *pattern, unsigned char current_pattern_length, unsigned int*result);
unsigned char search_pattern_in_FM_index_exSA_without_hash(char *pattern, unsigned char current_pattern_length, unsigned int*result);
long long int approximate_search_in_rankbitvector(char *pattern, unsigned int*result);
unsigned int get_SA_value_rankbitvector(unsigned int bwt_position);
unsigned char search_pattern_in_rankbitvector(char*pattern,unsigned char current_pattern_length, unsigned int*result);
unsigned int wt_rank_bitvector(unsigned char c, unsigned int position);
void create_rank_bitvectors(unsigned char*bwt);
long long int approximate_search_in_FM_index_entry(char *pattern, unsigned int*result);
void threshold_search_pattern_in_FM_index_entry(char *pattern, unsigned int *result_length, unsigned int *result, unsigned int pattern_length);
unsigned int extend_seed_in_FM_index_entry(unsigned char*pattern, unsigned int last, unsigned int*result);
unsigned char search_pattern_in_FM_index_entry(char *pattern, unsigned char current_pattern_length, unsigned int*result);
unsigned int get_SA_value_entry(unsigned int bwt_position);
unsigned char wt_access_entry(unsigned int entry_index, unsigned short int in_entry_index);
unsigned int wt_rank_entry(unsigned char c, unsigned int entry_index, unsigned int in_entry_index);
unsigned int wt_access_rank_entry(unsigned int entry_index, unsigned int in_entry_index);
void rebuild_FM_index_into_entries(unsigned int*suffix_array, unsigned char*bwt,unsigned char flag_withSA);
struct FMIndex_WT*build_FM_index_WT(unsigned int *suffix_array, unsigned char *bwt);
unsigned int find_end(unsigned int *suffix_array);
unsigned int *create_sample_SA(unsigned int *suffix_array);
unsigned int *create_count_table(unsigned char *s);
unsigned char get_alphabet_index(char *alphabet, unsigned char c);
unsigned int score(char a, char b);
unsigned char align(char *p1, char*p2);
unsigned int calculate(int i, int j, int k);
unsigned int get_max_array(unsigned char* array, unsigned int length);
unsigned int last_to_first_encoded(unsigned char c, unsigned int bwt_position, unsigned char *alphabet,unsigned int *count_table,unsigned char*bwt,unsigned int**occurence_table);
unsigned char search_pattern_in_FM_index_WT(char *pattern, unsigned char current_pattern_length, unsigned int*result);
long long int approximate_search_in_FM_index_WT(char *pattern, unsigned int*result);
unsigned int get_SA_value_WT(unsigned int bwt_position);
unsigned int last_to_first_WT(unsigned char c, unsigned int bwt_position, struct wavelet_tree *wtree_root,unsigned int *count_table);
void threshold_search_pattern_in_FM_index_WT(char *pattern, unsigned int *result_length, unsigned int *result,unsigned int pattern_length);
unsigned int extend_seed_in_FM_index_WT(unsigned char*pattern, unsigned int last, unsigned int*result);
void index_to_dna(unsigned int index, char *dna);
unsigned int dna_to_index(char *s);
unsigned char search_pattern_in_FM_index_entry_without_hash(char *pattern, unsigned char current_pattern_length, unsigned int*result);