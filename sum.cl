__kernel void sum(
	__global const char *input,
	__global unsigned int *output, 
	const int readLength, 
	__local unsigned int* start_end,
	__global const unsigned int* count_table_glob,
	__local unsigned int*count_table,
	__local char* alphabet_indexes,
	__local unsigned int* fm_index_entry,
	__local unsigned int* count_table_results,
	__local unsigned int* bitcounts,
	__global unsigned int*fm_index_input) 

{
	const int tid = get_local_id(0); //tid
	//const int globalID = BLOCK_SIZE; //i
	const int blockDim = get_local_size(0);//blockDim.x
	const int workgroupID = get_group_id(0); //blockIdx.x
	unsigned int entry_index;
	unsigned char in_entry_index = 0;
	unsigned int entry_index2;
	unsigned char in_entry_index2 = 0;
	
	int a;
	//const int gridSize = get_num_groups(0)*2*blockDim;
	//const int i = workgroupID*blockDim+tid;
	
	/*char warp_index = tid/32;
	char in_warp_index = tid%32;*/
	unsigned char iterator = readLength/2 - 1;
	int readIndex = workgroupID*readLength + tid + 32;

	barrier(CLK_LOCAL_MEM_FENCE);
	alphabet_indexes[tid] = input[readIndex];
	barrier(CLK_LOCAL_MEM_FENCE);
	alphabet_indexes[tid] = (alphabet_indexes[tid] == 'A') ? 0 : (alphabet_indexes[tid] == 'C') ? 1 : (alphabet_indexes[tid] == 'G') ? 2 : 3;



	//wait to load count_table
	barrier(CLK_LOCAL_MEM_FENCE);
	count_table_results[tid] = count_table_glob[alphabet_indexes[tid]];
	barrier(CLK_LOCAL_MEM_FENCE);
	if (tid == iterator)
		count_table_results[blockDim+tid] = count_table_glob[alphabet_indexes[tid] + 1] - 1;
	else
		count_table_results[blockDim+tid] = count_table_glob[alphabet_indexes[tid]];

	//wait to calc counts
	barrier(CLK_LOCAL_MEM_FENCE);
	
	while(iterator>0 && count_table_results[iterator] < count_table_results[iterator+blockDim])
	{

		entry_index = (count_table_results[iterator]>>8) * 32;
		entry_index2 = (count_table_results[iterator+blockDim]>>8) * 32;
		if (tid<8){
			in_entry_index = count_table_results[iterator] & 255;
			in_entry_index2 = count_table_results[iterator+blockDim] & 255;
		}


		iterator--;
		barrier(CLK_LOCAL_MEM_FENCE);
		//loading fm index entry
		fm_index_entry[tid] = fm_index_input[ entry_index + tid];
		barrier(CLK_LOCAL_MEM_FENCE);

		////pocitam dalsi rozsah, teraz spodna hranica
		//max prvych 8 threadov z 32
		if (tid<8)
		{
			if (tid<(in_entry_index >> 5))
			{	
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}
			
			else if (tid == (in_entry_index >> 5))
			{	
				if (in_entry_index & 31)
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - (in_entry_index & 31))); 
				else
					bitcounts[tid] = 0;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (tid == 0)
		{	
			for (a=1;a<=(in_entry_index >> 5);a++)
				bitcounts[tid] += bitcounts[a];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (tid == 0 && (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3))
				bitcounts[tid] = in_entry_index - bitcounts[tid];

		barrier(CLK_LOCAL_MEM_FENCE);
		//max 8 threadov, 9-17 pozicia
		//spocitat 'vlavo"
		if (tid >  8 && (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3))
		{
			if (tid<((bitcounts[0]/32) + 9))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}
		
			else if (tid == (9+bitcounts[0]/32))
			{	
				if (bitcounts[0]%32)
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - bitcounts[0]%32)); 
				else
					bitcounts[tid] = 0;
			}
		}
		//spocitat "vpravo"
		if (tid >  8 && (alphabet_indexes[iterator] == 1 || alphabet_indexes[iterator] == 2))
		{
			a = fm_index_entry[8];

			if (tid<(9+bitcounts[0]/32))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid + a]);
			}
		
			else if (tid == (9+bitcounts[0]/32))
			{	
				if (bitcounts[0]%32){
					bitcounts[tid] = popcount((fm_index_entry[tid + a]) >> (32 - bitcounts[0]%32)); 

				}
				else 
					bitcounts[tid] = 0;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (tid == 9)
		{
			if (alphabet_indexes[iterator] == 3 || alphabet_indexes[iterator] == 1)
			{
				for (a=1;a<=bitcounts[0]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
			}	
			else
			{
				for (a=1;a<=bitcounts[0]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
				bitcounts[tid] = bitcounts[0] - bitcounts[tid];

			}
		}

		//add occ_counter in fm index entry
		if (tid == iterator)
		{
			if (alphabet_indexes[iterator] == 0)
				count_table_results[iterator] += fm_index_entry[20];
			else if (alphabet_indexes[iterator] == 1)
				count_table_results[iterator] += fm_index_entry[23];
			else if (alphabet_indexes[iterator] == 2)
				count_table_results[iterator] += fm_index_entry[22];
			else if (alphabet_indexes[iterator] == 3)
				count_table_results[iterator] += fm_index_entry[21];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (tid == 9){
			//printf("pred %d a scita sa %d\n",count_table_results[iterator], bitcounts[tid]);
			count_table_results[iterator] += bitcounts[tid]; //fm_index_occ is already loaded
		}

		////pocitam dalsi rozsah, teraz horna hranicass
		//loading fm index entry


		barrier(CLK_LOCAL_MEM_FENCE);
		fm_index_entry[tid] = fm_index_input[ entry_index2 + tid];
		barrier(CLK_LOCAL_MEM_FENCE);
		
		if (tid<8)
		{
			if (tid<in_entry_index2/32)
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			
			else if (tid == in_entry_index2/32)
			{
				if (in_entry_index2%32)
					bitcounts[tid] = popcount((fm_index_entry[tid]) >> (32 - in_entry_index2%32)); 
				else
					bitcounts[tid] = 0;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (tid == 0)
		{
			for (a=1;a<=in_entry_index2/32;a++)
				bitcounts[tid] += bitcounts[a];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (tid ==0 && (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3))
				bitcounts[tid] = in_entry_index2 - bitcounts[0];
			
		barrier(CLK_LOCAL_MEM_FENCE);
		//max 8 threadov, 9-17 pozicia
		//spocitat 'vlavo"
		if (tid >  8 && (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3))
		{
			if (tid<(9+bitcounts[0]/32))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}	
			else if (tid == (9+bitcounts[0]/32))
			{	
				if (bitcounts[0]%32){
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - bitcounts[0]%32)); 
				}
				else
					bitcounts[tid] = 0;
			}
		}
		//spocitat "vpravo"
		if (tid >  8 && (alphabet_indexes[iterator] == 1 || alphabet_indexes[iterator] == 2))
		{
			a = fm_index_entry[8];

			if (tid<(9+bitcounts[0]/32))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid + a]);
			}
			else if (tid == (9+bitcounts[0]/32))
			{	
				if (bitcounts[0]%32){
					bitcounts[tid] = popcount((fm_index_entry[tid + a]) >> (32 - bitcounts[0]%32)); 
				}
				else 
					bitcounts[tid] = 0;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (tid == 9)
		{
			if (alphabet_indexes[iterator] == 3 || alphabet_indexes[iterator] == 1)
			{
				for (a=1;a<=bitcounts[0]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
			}	
			else
			{
				for (a=1;a<=bitcounts[0]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
				bitcounts[tid] = bitcounts[0] - bitcounts[tid];
			}
		}
		//add occ_counter in fm index entry
		if (tid == 0)
		{
			if (alphabet_indexes[iterator] == 0)
				count_table_results[iterator+blockDim] += fm_index_entry[20];
			else if (alphabet_indexes[iterator] == 1)
				count_table_results[iterator+blockDim] += fm_index_entry[23];
			else if (alphabet_indexes[iterator] == 2)
				count_table_results[iterator+blockDim] += fm_index_entry[22];
			else if (alphabet_indexes[iterator] == 3)
				count_table_results[iterator+blockDim] += fm_index_entry[21];
			
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (tid == 9)
		{
			count_table_results[iterator+blockDim] += bitcounts[tid]; //fm_index_occ is already loaded
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	//KONIEC VYHLADANIA CASTI

	
	if (tid ==0 ){
		output[4*workgroupID] = count_table_results[iterator];
		output[4*workgroupID+1] = count_table_results[iterator+blockDim];
	}



	iterator = readLength/2 - 1;
	readIndex = readIndex - 32;
	alphabet_indexes[tid] = input[readIndex];
	alphabet_indexes[tid] = (alphabet_indexes[tid] == 'A') ? 0 : (alphabet_indexes[tid] == 'C') ? 1 : (alphabet_indexes[tid] == 'G') ? 2 : 3;

	count_table_results[tid] = count_table_glob[alphabet_indexes[tid]];
	barrier(CLK_LOCAL_MEM_FENCE);
	if (tid == iterator)
		count_table_results[blockDim+tid] = count_table_glob[alphabet_indexes[tid] + 1] - 1;
	else
		count_table_results[blockDim+tid] = count_table_glob[alphabet_indexes[tid]];

	//wait to calc counts
	barrier(CLK_LOCAL_MEM_FENCE);
	
	
	while(iterator>0 && count_table_results[iterator] < count_table_results[iterator+blockDim])
	{

		entry_index = (count_table_results[iterator]>>8) * 32;
		entry_index2 = (count_table_results[iterator+blockDim]>>8) * 32;
		if (tid<8){
			in_entry_index = count_table_results[iterator] & 255;
			in_entry_index2 = count_table_results[iterator+blockDim] & 255;
		}

		barrier(CLK_LOCAL_MEM_FENCE);

		iterator--;
		barrier(CLK_LOCAL_MEM_FENCE);
		//loading fm index entry
		fm_index_entry[tid] = fm_index_input[ entry_index + tid];
		barrier(CLK_LOCAL_MEM_FENCE);

		////pocitam dalsi rozsah, teraz spodna hranica
		//max prvych 8 threadov z 32
		if (tid<8)
		{
			if (tid<(in_entry_index >> 5))
			{	
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}
			
			else if (tid == (in_entry_index >> 5))
			{	
				if ((in_entry_index & 31))
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - (in_entry_index & 31))); 
				else
					bitcounts[tid] = 0;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (tid == 0)
		{	
			for (a=1;a<=(in_entry_index >> 5);a++)
				bitcounts[tid] += bitcounts[a];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3)
		{
			if (tid == 0)
			{	
				bitcounts[tid] = in_entry_index - bitcounts[tid];
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//max 8 threadov, 9-17 pozicia
		//spocitat 'vlavo"
		if (tid >  8 && (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3))
		{
			if (tid<((bitcounts[0]/32) + 9))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}
		
			else if (tid == (9+bitcounts[0]/32))
			{	
				if (bitcounts[0]%32)
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - bitcounts[0]%32)); 
				else
					bitcounts[tid] = 0;
			}
		}
		//spocitat "vpravo"
		if (tid >  8 && (alphabet_indexes[iterator] == 1 || alphabet_indexes[iterator] == 2))
		{
			a = fm_index_entry[8];

			if (tid<(9+bitcounts[0]/32))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid + a]);
			}
		
			else if (tid == (9+bitcounts[0]/32))
			{	
				if (bitcounts[0]%32){
					bitcounts[tid] = popcount((fm_index_entry[tid + a]) >> (32 - bitcounts[0]%32)); 

				}
				else 
					bitcounts[tid] = 0;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (tid == 9)
		{
			if (alphabet_indexes[iterator] == 3 || alphabet_indexes[iterator] == 1)
			{
				for (a=1;a<=bitcounts[0]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
			}	
			else
			{
				for (a=1;a<=bitcounts[0]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
				bitcounts[tid] = bitcounts[0] - bitcounts[tid];

			}
		}

		//add occ_counter in fm index entry
		else if (tid == 0)
		{
			if (alphabet_indexes[iterator] == 0)
				count_table_results[iterator] += fm_index_entry[20];
			else if (alphabet_indexes[iterator] == 1)
				count_table_results[iterator] += fm_index_entry[23];
			else if (alphabet_indexes[iterator] == 2)
				count_table_results[iterator] += fm_index_entry[22];
			else if (alphabet_indexes[iterator] == 3)
				count_table_results[iterator] += fm_index_entry[21];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (tid == 9){
			//printf("pred %d a scita sa %d\n",count_table_results[iterator], bitcounts[tid]);
			count_table_results[iterator] += bitcounts[tid]; //fm_index_occ is already loaded
		}

		////pocitam dalsi rozsah, teraz horna hranicass
		//loading fm index entry


		barrier(CLK_LOCAL_MEM_FENCE);
		fm_index_entry[tid] = fm_index_input[ entry_index2 + tid];
		barrier(CLK_LOCAL_MEM_FENCE);
		
		if (tid<in_entry_index2/32)
		{
			bitcounts[tid] = popcount(fm_index_entry[tid]);
		}
		
		else if (tid == in_entry_index2/32)
		{
			if (in_entry_index2%32)
				bitcounts[tid] = popcount((fm_index_entry[tid]) >> (32 - in_entry_index2%32)); 
			else
				bitcounts[tid] = 0;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (tid == 0)
		{
			for (a=1;a<=in_entry_index2/32;a++)
				bitcounts[tid] += bitcounts[a];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3)
		{
			if (tid == 0)
			{
				bitcounts[tid] = in_entry_index2 - bitcounts[0];
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//max 8 threadov, 9-17 pozicia
		//spocitat 'vlavo"
		if (tid >  8 && (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3))
		{
			if (tid<(9+bitcounts[0]/32))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}	
			else if (tid == (9+bitcounts[0]/32))
			{	
				if (bitcounts[0]%32){
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - bitcounts[0]%32)); 
				}
				else
					bitcounts[tid] = 0;
			}
		}
		//spocitat "vpravo"
		if (tid >  8 && (alphabet_indexes[iterator] == 1 || alphabet_indexes[iterator] == 2))
		{
			a = fm_index_entry[8];

			if (tid<(9+bitcounts[0]/32))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid + a]);
			}
			else if (tid == (9+bitcounts[0]/32))
			{	
				if (bitcounts[0]%32){
					bitcounts[tid] = popcount((fm_index_entry[tid + a]) >> (32 - bitcounts[0]%32)); 
				}
				else 
					bitcounts[tid] = 0;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (tid == 9)
		{
			if (alphabet_indexes[iterator] == 3 || alphabet_indexes[iterator] == 1)
			{
				for (a=1;a<=bitcounts[0]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
			}	
			else
			{
				for (a=1;a<=bitcounts[0]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
				bitcounts[tid] = bitcounts[0] - bitcounts[tid];
			}
		}
		//add occ_counter in fm index entry
		if (tid == 0)
		{
			if (alphabet_indexes[iterator] == 0)
				count_table_results[iterator+blockDim] += fm_index_entry[20];
			else if (alphabet_indexes[iterator] == 1)
				count_table_results[iterator+blockDim] += fm_index_entry[23];
			else if (alphabet_indexes[iterator] == 2)
				count_table_results[iterator+blockDim] += fm_index_entry[22];
			else if (alphabet_indexes[iterator] == 3)
				count_table_results[iterator+blockDim] += fm_index_entry[21];
			
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (tid == 9)
		{
			count_table_results[iterator+blockDim] += bitcounts[tid]; //fm_index_occ is already loaded
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	//KONIEC VYHLADANIA CASTI

	if (tid ==0 ){
		output[4*workgroupID+2] = count_table_results[iterator];
		output[4*workgroupID+3] = count_table_results[iterator+blockDim];
	}
}
	
