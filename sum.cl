__kernel void sum(
	__global const char *input,
	__global unsigned int *output, 
	const int readLength, 
	__local unsigned int* start_end,
	__global const unsigned int* count_table_glob,
	__local unsigned int*count_table,
	__local char* alphabet_indexes,
	__local int* indexes,
	__local unsigned int* fm_index_entry,
	__local unsigned int* count_table_results,
	__local unsigned int* bitcounts,
	__global unsigned int*fm_index_input) 

{
	const int tid = get_local_id(0); //tid
	//const int globalID = BLOCK_SIZE; //i
	const int blockDim = get_local_size(0);//blockDim.x
	const int workgroupID = get_group_id(0); //blockIdx.x
	int i;
	unsigned int temp;
	

	int a;
	const int gridSize = get_num_groups(0)*2*blockDim;
	//const int i = workgroupID*blockDim+tid;
	
	
	/*char warp_index = tid/32;
	char in_warp_index = tid%32;*/
	unsigned char iterator = readLength/2 - 1;

	int readIndex = workgroupID*readLength + tid + 32;

	 /*If n is a power of 2, (i/n) is equivalent to (i ≫ log2(n)) 
	 and (i % n) is equivalent to (i & (n-1)).*/


	//loading read ...
	//for (a=0;a<readLength/32;a++)
	
	
	alphabet_indexes[tid] = (input[readIndex] == 'A') ? 0 : (input[readIndex] == 'C') ? 1 : (input[readIndex] == 'G') ? 2 : 3;
	//printf("%d vlakno, tid %d, blockdim %d, workID %d, read %d, iter %d po %d, c %c, ind%d\n",i,tid,blockDim,workgroupID,readIndex,iterator, warp_index*32, input[readIndex], alphabet_indexes[tid]);
	
	/*readIndex = readIndex - 32;
	alphabet_indexes[tid] = (input[readIndex] == 'A') ? 0 : (input[readIndex] == 'C') ? 1 : (input[readIndex] == 'G') ? 2 : 3;
	*/
	

	//loading count_table global? IDK IF NEEDED WHEN CONSTANT
	if (tid<5)
	{
		count_table[tid] = count_table_glob[tid];
	}

	//wait to load count_table
	barrier(CLK_LOCAL_MEM_FENCE);
	count_table_results[tid] = count_table[alphabet_indexes[tid]];
	barrier(CLK_LOCAL_MEM_FENCE);
	if (tid == iterator)
		count_table_results[blockDim+tid] = count_table[alphabet_indexes[tid] + 1] - 1;
	else
		count_table_results[blockDim+tid] = count_table[alphabet_indexes[tid]];

	//wait to calc counts
	barrier(CLK_LOCAL_MEM_FENCE);
	
	while(iterator>0 && count_table_results[iterator] < count_table_results[iterator+blockDim])
	{
		if (tid<4)
		{
			if (tid == 0)
			{
			 indexes[0] = (count_table_results[iterator]>>8) * 32;	//entry index
			}
			else if (tid == 1)
			{
			 indexes[1] = count_table_results[iterator] & 255;	//in_entry_index
			}
			else if (tid == 2)
			{
			 indexes[2] = (count_table_results[blockDim+iterator]>>8) * 32;	//entry index
			}
			else
			{	
			 indexes[3] = count_table_results[blockDim+iterator]&255;	//in_entry_index
			}
		}

		barrier(CLK_LOCAL_MEM_FENCE);
		if (tid>=4 && tid<8)
		{
			if (tid == 4)
			{
			 indexes[4] = indexes[1] >> 5;	//entry index
			}
			else if (tid == 5)
			{
			 indexes[5] = indexes[1] & 31;	//in_entry_index
			}
			else if (tid == 6)
			{
			 indexes[6] = indexes[3] >> 5;	//entry_index
			}
			else if (tid == 7)
			{	
			 indexes[7] = indexes[3] & 31;	//in_entry_index
			}
		}

		iterator--;
		barrier(CLK_LOCAL_MEM_FENCE);
		//loading fm index entry
		fm_index_entry[tid] = fm_index_input[ indexes[0] + tid];
		barrier(CLK_LOCAL_MEM_FENCE);

		////pocitam dalsi rozsah, teraz spodna hranica
		//max prvych 8 threadov z 32
		if (tid<8)
		{
			if (tid<indexes[4])
			{	
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}
			
			else if (tid == indexes[4])
			{	
				if (indexes[5])
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - indexes[5])); 
				else
					bitcounts[tid] = 0;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (tid == 0)
		{	
			for (a=1;a<=indexes[4];a++)
				bitcounts[tid] += bitcounts[a];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3)
		{
			if (tid == 0)
			{	
				bitcounts[tid] = indexes[1] - bitcounts[tid];
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
		fm_index_entry[tid] = fm_index_input[ indexes[2] + tid];
		barrier(CLK_LOCAL_MEM_FENCE);
		
		if (tid<indexes[3]/32)
		{
			bitcounts[tid] = popcount(fm_index_entry[tid]);
		}
		
		else if (tid == indexes[3]/32)
		{
			if (indexes[3]%32)
				bitcounts[tid] = popcount((fm_index_entry[tid]) >> (32 - indexes[3]%32)); 
			else
				bitcounts[tid] = 0;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (tid == 0)
		{
			for (a=1;a<=indexes[3]/32;a++)
				bitcounts[tid] += bitcounts[a];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3)
		{
			if (tid == 0)
			{
				bitcounts[tid] = indexes[3] - bitcounts[0];
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
		output[4*workgroupID] = count_table_results[iterator];
		output[4*workgroupID+1] = count_table_results[iterator+blockDim];
	}



	iterator = readLength/2 - 1;
	readIndex = readIndex - 32;
	alphabet_indexes[tid] = (input[readIndex] == 'A') ? 0 : (input[readIndex] == 'C') ? 1 : (input[readIndex] == 'G') ? 2 : 3;

	count_table_results[tid] = count_table[alphabet_indexes[tid]];
	barrier(CLK_LOCAL_MEM_FENCE);
	if (tid == iterator)
		count_table_results[blockDim+tid] = count_table[alphabet_indexes[tid] + 1] - 1;
	else
		count_table_results[blockDim+tid] = count_table[alphabet_indexes[tid]];

	//wait to calc counts
	barrier(CLK_LOCAL_MEM_FENCE);
	
	

	while(iterator>0 && count_table_results[iterator] < count_table_results[iterator+blockDim])
	{

		if (tid<4)
		{
			if (tid == 0)
			{
			 indexes[0] = (count_table_results[iterator]>>8) * 32;	//entry index
			}
			else if (tid == 1)
			{
			 indexes[1] = count_table_results[iterator] & 255;	//in_entry_index
			}
			else if (tid == 2)
			{
			 indexes[2] = (count_table_results[blockDim+iterator]>>8) * 32;	//entry index
			}
			else
			{	
			 indexes[3] = count_table_results[blockDim+iterator]&255;	//in_entry_index
			}
		}

		barrier(CLK_LOCAL_MEM_FENCE);
		if (tid>=4 && tid<8)
		{
			if (tid == 4)
			{
			 indexes[4] = indexes[1] >> 5;	//entry index
			}
			else if (tid == 5)
			{
			 indexes[5] = indexes[1] & 31;	//in_entry_index
			}
			else if (tid == 6)
			{
			 indexes[6] = indexes[3] >> 5;	//entry_index
			}
			else if (tid == 7)
			{	
			 indexes[7] = indexes[3] & 31;	//in_entry_index
			}
		}

		iterator--;
		barrier(CLK_LOCAL_MEM_FENCE);
		//loading fm index entry
		fm_index_entry[tid] = fm_index_input[ indexes[0] + tid];
		barrier(CLK_LOCAL_MEM_FENCE);

		////pocitam dalsi rozsah, teraz spodna hranica
		//max prvych 8 threadov z 32
		if (tid<8)
		{
			if (tid<indexes[4])
			{	
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}
			
			else if (tid == indexes[4])
			{	
				if (indexes[5])
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - indexes[5])); 
				else
					bitcounts[tid] = 0;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (tid == 0)
		{	
			for (a=1;a<=indexes[4];a++)
				bitcounts[tid] += bitcounts[a];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3)
		{
			if (tid == 0)
			{	
				bitcounts[tid] = indexes[1] - bitcounts[tid];
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
		fm_index_entry[tid] = fm_index_input[ indexes[2] + tid];
		barrier(CLK_LOCAL_MEM_FENCE);
		
		if (tid<indexes[3]/32)
		{
			bitcounts[tid] = popcount(fm_index_entry[tid]);
		}
		
		else if (tid == indexes[3]/32)
		{
			if (indexes[3]%32)
				bitcounts[tid] = popcount((fm_index_entry[tid]) >> (32 - indexes[3]%32)); 
			else
				bitcounts[tid] = 0;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (tid == 0)
		{
			for (a=1;a<=indexes[3]/32;a++)
				bitcounts[tid] += bitcounts[a];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3)
		{
			if (tid == 0)
			{
				bitcounts[tid] = indexes[3] - bitcounts[0];
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
	

	/*temp = count_table_results[iterator];
	while( temp < count_table_results[iterator+blockDim])
	{

		count_table_results[iterator] = temp;

		i = 0;	//na pocitadlo
		while (count_table_results[iterator]%32 != 0)
		{

			if (tid == 1)
				indexes[0] = (count_table_results[iterator]/256)*32;
			else if (tid == 2)
				indexes[1] = count_table_results[iterator]%256;


			barrier(CLK_LOCAL_MEM_FENCE);
			
			fm_index_entry[tid] = fm_index_input[ indexes[0] + tid];
			
			barrier(CLK_LOCAL_MEM_FENCE);

		


			//wt_access
			//count 0 or 1 up to indexes[2*warp_index+1]
			//find 

			if (tid<8)
			{
				if (tid<indexes[1]/32)
				{	
					bitcounts[tid] = popcount(fm_index_entry[tid]);
				}
				
				else if (tid == indexes[1]/32)
				{	
					if (indexes[1]%32)
						bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - (indexes[1]%32))); 
					else
						bitcounts[tid] = 0;
				}
			}
			else if (tid == 8)
				bitcounts[tid] = 31 - indexes[1]%32;
			else if (tid == 9)
				bitcounts[tid] = indexes[1]/32;

			barrier(CLK_LOCAL_MEM_FENCE);
			//counting together
			if (tid == 0)
			{
				for (a=1;a<=indexes[1]/32;a++)
					bitcounts[tid] += bitcounts[a];
			}

			else if (tid == 8) // v 8 ulozim aky je bit
			{
				bitcounts[tid] = (fm_index_entry[bitcounts[tid+1]] >> bitcounts[tid]) & 1;
			}
			barrier(CLK_LOCAL_MEM_FENCE);

			//v 0.tej pozicii je suma 1
			//v 1.pozicii je uchovana suma 0
			//v 2. pozicii je suma %32
			//v 3.pozicii je suma/32
			if (tid == 1)
			{ 
				if (bitcounts[8])
				{
					bitcounts[tid] = bitcounts[0];
				}
				else
					bitcounts[tid] = indexes[1] - bitcounts[0];
			}

			if (tid == 0 && workgroupID == 0)
				printf("sum1 %d, sum0 %d, zvy %d, pod %d, c %d\n",bitcounts[0],bitcounts[1], bitcounts[2], bitcounts[3], bitcounts[8]);

			barrier(CLK_LOCAL_MEM_FENCE);
			if (tid == 2)
				bitcounts[tid] = 31  - bitcounts[1]%32;
			else if (tid == 3)
				bitcounts[tid] = bitcounts[1]/32;
			
			barrier(CLK_LOCAL_MEM_FENCE);
			if (tid == 9 )
			{
				if (bitcounts[8])
				{
					if (((fm_index_entry[bitcounts[3] + fm_index_entry[8] + 9]) >> bitcounts[2]) & 1)
						alphabet_indexes[iterator] = 1;
					else
						alphabet_indexes[iterator] = 2;
				}
				else
				{
					if (((fm_index_entry[bitcounts[3] + 9]) >> bitcounts[2]) & 1)
						alphabet_indexes[iterator] = 3;
					else
						alphabet_indexes[iterator] = 0;
				}
			}		
			barrier(CLK_LOCAL_MEM_FENCE);

			//koniec wt_access
			//1 z 32
			if (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3)
			{
				if (tid == 0)
				{
					bitcounts[tid] = indexes[1] - bitcounts[0];
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
					count_table_results[iterator] = count_table[0] + fm_index_entry[20];
				else if (alphabet_indexes[iterator] == 1)
					count_table_results[iterator] = count_table[1] + fm_index_entry[23];
				else if (alphabet_indexes[iterator] == 2)
					count_table_results[iterator] = count_table[2] + fm_index_entry[22];
				else if (alphabet_indexes[iterator] == 3)
					count_table_results[iterator] = count_table[3] + fm_index_entry[21];
			}
			barrier(CLK_LOCAL_MEM_FENCE);
			if (tid == 9)
			{
				count_table_results[iterator] += bitcounts[tid]; //fm_index_occ is already loaded
			}
			i++;
			barrier(CLK_LOCAL_MEM_FENCE);
		}

		
		if (tid == 0 && workgroupID ==0)
			printf("sa result je %d, count je %d\n",count_table_results[iterator],i);
		temp++;
	}
	



	if (tid ==0 ){
		printf("vraciam %d\n",workgroupID);

		output[workgroupID] = 33;
	}
}
	*/