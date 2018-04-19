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
	
	
	char warp_index = tid/32;
	char in_warp_index = tid%32;
	unsigned char iterator = readLength/2 - 1 + warp_index*32;

	int readIndex = readLength * warp_index + workgroupID*8*readLength + in_warp_index + 32;

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


	if (tid == iterator)
	{
		count_table_results[tid] = count_table[alphabet_indexes[tid]];
		count_table_results[blockDim+tid] = count_table[alphabet_indexes[tid] + 1] - 1;

	}
	else
	{
		count_table_results[tid] = count_table[alphabet_indexes[tid]];
		count_table_results[blockDim+tid] = count_table[alphabet_indexes[tid]];
	}


	
	//wait to calc counts
	barrier(CLK_LOCAL_MEM_FENCE);
	

	//iterator--;
	
	while(iterator>warp_index*32 && count_table_results[iterator] < count_table_results[iterator+blockDim])
	{


		if (in_warp_index == 0)
		{
		 indexes[warp_index*4] = (count_table_results[iterator]/256) * 32;	//entry index
		}
		else if (in_warp_index == 1)
		{
		 indexes[warp_index*4 + 1] = count_table_results[iterator]%256;	//in_entry_index
		}
		else if (in_warp_index == 2)
		{
		 indexes[warp_index*4 + 2] = (count_table_results[256+iterator]/256) * 32;	//entry index
		}
		else if (in_warp_index == 3)
		{	
		 indexes[warp_index*4 + 3] = count_table_results[256+iterator]%256;	//in_entry_index
		}


		iterator--;
		
		barrier(CLK_LOCAL_MEM_FENCE);

		//loading fm index entry
		fm_index_entry[tid] = fm_index_input[ indexes[4*warp_index] + in_warp_index];


		barrier(CLK_LOCAL_MEM_FENCE);


		////pocitam dalsi rozsah, teraz spodna hranica
		//max prvych 8 threadov z 32
		if (in_warp_index<8)
		{
			if (in_warp_index<indexes[warp_index*4 + 1]/32)
			{	
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}
			
			else if (in_warp_index == indexes[warp_index*4 + 1]/32)
			{	
				if (indexes[warp_index*4 + 1]%32)
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - (indexes[warp_index*4 + 1]%32))); 
				else
					bitcounts[tid] = 0;
				
				//printf("je %d hodnot, pocita akurat %d vlakno so zvyskm, vrati %d, shift %d,  vnutri %d, po shifte %lu\n",indexes[warp_index*4 + 1], in_warp_index, bitcounts[tid], (32 - (indexes[warp_index*4 + 1]%32)),indexes[warp_index*4 + 1]%32, fm_index_entry[tid]>>32);
			}
		}

		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (in_warp_index == 0)
		{	
			for (a=1;a<=indexes[warp_index*4 + 1]/32;a++)
				bitcounts[tid] += bitcounts[tid + a];
		}


		barrier(CLK_LOCAL_MEM_FENCE);

		//1 z 32
		if (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3)
		{
			if (in_warp_index == 0)
			{	
				bitcounts[tid] = indexes[warp_index*4 + 1] - bitcounts[warp_index*32];
			}
		}
		
		barrier(CLK_LOCAL_MEM_FENCE);
		
		//max 8 threadov, 9-17 pozicia
		//spocitat 'vlavo"
		if (in_warp_index >  8 && (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3))
		{
			if (in_warp_index<(9+bitcounts[warp_index*32]/32))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}
		
			else if (in_warp_index == (9+bitcounts[warp_index*32]/32))
			{	
				if (bitcounts[warp_index*32]%32)
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - bitcounts[warp_index*32]%32)); 
				else
					bitcounts[tid] = 0;
			}
		}

		//spocitat "vpravo"
		if (in_warp_index >  8 && (alphabet_indexes[iterator] == 1 || alphabet_indexes[iterator] == 2))
		{
			a = fm_index_entry[warp_index*32 + 8];

			if (in_warp_index<(9+bitcounts[warp_index*32]/32))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid + a]);
			}
		
			else if (in_warp_index == (9+bitcounts[warp_index*32]/32))
			{	
				if (bitcounts[warp_index*32]%32){
					bitcounts[tid] = popcount((fm_index_entry[tid + a]) >> (32 - bitcounts[warp_index*32]%32)); 

				}
				else 
					bitcounts[tid] = 0;
			}
		}


		barrier(CLK_LOCAL_MEM_FENCE);

		if (in_warp_index == 9)
		{
			if (alphabet_indexes[iterator] == 3 || alphabet_indexes[iterator] == 1)
			{
				for (a=1;a<=bitcounts[warp_index*32]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
			}	
			else
			{
				for (a=1;a<=bitcounts[warp_index*32]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
				bitcounts[tid] = bitcounts[warp_index*32] - bitcounts[tid];

			}
		}

		

		//add occ_counter in fm index entry
		if (in_warp_index == 0)
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

		if (in_warp_index == 9){
			//printf("pred %d a scita sa %d\n",count_table_results[iterator], bitcounts[tid]);
			count_table_results[iterator] += bitcounts[tid]; //fm_index_occ is already loaded
		}




		////pocitam dalsi rozsah, teraz horna hranicass
		//loading fm index entry
		barrier(CLK_LOCAL_MEM_FENCE);
		fm_index_entry[tid] = fm_index_input[ indexes[4*warp_index + 2] + in_warp_index];
		barrier(CLK_LOCAL_MEM_FENCE);
		

		if (in_warp_index<indexes[warp_index*4 + 3]/32)
		{
			bitcounts[tid] = popcount(fm_index_entry[tid]);
		}
		
		else if (in_warp_index == indexes[warp_index*4 + 3]/32)
		{
			if (indexes[warp_index*4 + 3]%32)
				bitcounts[tid] = popcount((fm_index_entry[tid]) >> (32 - indexes[warp_index*4 + 3]%32)); 
			else
				bitcounts[tid] = 0;
		}



		barrier(CLK_LOCAL_MEM_FENCE);
		//1 z 32
		if (in_warp_index == 0)
		{
			for (a=1;a<=indexes[warp_index*4 + 3]/32;a++)
				bitcounts[tid] += bitcounts[tid + a];
		}

		barrier(CLK_LOCAL_MEM_FENCE);

		//1 z 32
		if (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3)
		{
			if (in_warp_index == 0)
			{
				bitcounts[tid] = indexes[warp_index*4 + 3] - bitcounts[warp_index*32];
			}
		}
		

		barrier(CLK_LOCAL_MEM_FENCE);


		//max 8 threadov, 9-17 pozicia
		//spocitat 'vlavo"
		if (in_warp_index >  8 && (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3))
		{

			if (in_warp_index<(9+bitcounts[warp_index*32]/32))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid]);
			}
		
			else if (in_warp_index == (9+bitcounts[warp_index*32]/32))
			{	
				if (bitcounts[warp_index*32]%32){
					bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - bitcounts[warp_index*32]%32)); 
				}
				else
					bitcounts[tid] = 0;
			}
		}

		//spocitat "vpravo"
		if (in_warp_index >  8 && (alphabet_indexes[iterator] == 1 || alphabet_indexes[iterator] == 2))
		{

			a = fm_index_entry[warp_index*32 + 8];

			if (in_warp_index<(9+bitcounts[warp_index*32]/32))
			{
				bitcounts[tid] = popcount(fm_index_entry[tid + a]);
			}
		
			else if (in_warp_index == (9+bitcounts[warp_index*32]/32))
			{	
				if (bitcounts[warp_index*32]%32){
					bitcounts[tid] = popcount((fm_index_entry[tid + a]) >> (32 - bitcounts[warp_index*32]%32)); 
				}
				else 
					bitcounts[tid] = 0;
			}

			
		}

		barrier(CLK_LOCAL_MEM_FENCE);

		if (in_warp_index == 9)
		{
			if (alphabet_indexes[iterator] == 3 || alphabet_indexes[iterator] == 1)
			{
				for (a=1;a<=bitcounts[warp_index*32]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
			}	
			else
			{
				for (a=1;a<=bitcounts[warp_index*32]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
				bitcounts[tid] = bitcounts[warp_index*32] - bitcounts[tid];
			}
		}


		//add occ_counter in fm index entry
		if (in_warp_index == 0)
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

		if (in_warp_index == 9)
		{
			count_table_results[iterator+blockDim] += bitcounts[tid]; //fm_index_occ is already loaded
		}


		barrier(CLK_LOCAL_MEM_FENCE);


		/*if (tid==0 && workgroupID == 0)
	 		printf("iterator e %d a mal skoncit na %d, county %d %d\n",iterator, warp_index*32, count_table_results[iterator], count_table_results[iterator+blockDim]);
*/
	}

	temp = count_table_results[iterator];
	while( temp < count_table_results[iterator+blockDim])
	{

		count_table_results[iterator] = temp;

		i = 0;	//na pocitadlo
		while (count_table_results[iterator]%32 != 0)
		{

			if (in_warp_index == 1)
				indexes[warp_index*2] = (count_table_results[iterator]/256)*32;
			else if (in_warp_index == 2)
				indexes[warp_index*2 + 1] = count_table_results[iterator]%256;


			barrier(CLK_LOCAL_MEM_FENCE);
			
			fm_index_entry[tid] = fm_index_input[ indexes[2*warp_index] + in_warp_index];
			
			barrier(CLK_LOCAL_MEM_FENCE);

		


			//wt_access
			//count 0 or 1 up to indexes[2*warp_index+1]
			//find 

			if (in_warp_index<8)
			{
				if (in_warp_index<indexes[warp_index*2 + 1]/32)
				{	
					bitcounts[tid] = popcount(fm_index_entry[tid]);
				}
				
				else if (in_warp_index == indexes[warp_index*2 + 1]/32)
				{	
					if (indexes[warp_index*2 + 1]%32)
						bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - (indexes[warp_index*2 + 1]%32))); 
					else
						bitcounts[tid] = 0;
				}
			}

			else if (in_warp_index == 8)
				bitcounts[tid] = 31 - indexes[2*warp_index+1]%32;
			else if (in_warp_index == 9)
				bitcounts[tid] = indexes[2*warp_index+1]/32;


			barrier(CLK_LOCAL_MEM_FENCE);
			//counting together
			if (in_warp_index == 0)
			{
				for (a=1;a<=indexes[warp_index*2 + 1]/32;a++)
					bitcounts[tid] += bitcounts[tid + a];
			}

			else if (in_warp_index == 8) // v 8 ulozim aky je bit
			{
				bitcounts[tid] = (fm_index_entry[bitcounts[tid+1]] >> bitcounts[tid]) & 1;
			}
			barrier(CLK_LOCAL_MEM_FENCE);

			//v 0.tej pozicii je suma 1
			//v 1.pozicii je uchovana suma 0
			//v 2. pozicii je suma %32
			//v 3.pozicii je suma/32
			if (in_warp_index == 1)
			{ 
				if (bitcounts[warp_index*32 + 8])
				{
					bitcounts[tid] = bitcounts[warp_index*32];
				}
				else
					bitcounts[tid] = indexes[warp_index*2 + 1] - bitcounts[warp_index*32];
			}


			if (tid == 0 && workgroupID == 0)
				printf("sum1 %d, sum0 %d, zvy %d, pod %d, c %d\n",bitcounts[0],bitcounts[1], bitcounts[2], bitcounts[3], bitcounts[8]);

			barrier(CLK_LOCAL_MEM_FENCE);

			if (in_warp_index == 2)
				bitcounts[tid] = 31  - bitcounts[warp_index*32+1]%32;
			else if (in_warp_index == 3)
				bitcounts[tid] = bitcounts[warp_index*32+1]/32;
			
			barrier(CLK_LOCAL_MEM_FENCE);

			if (in_warp_index == 9 )
			{
				if (bitcounts[warp_index*32+8])
				{
					if (((fm_index_entry[bitcounts[warp_index*32 + 3] + fm_index_entry[warp_index*32 + 8] + 9]) >> bitcounts[warp_index*32 + 2]) & 1)
						alphabet_indexes[iterator] = 1;
					else
						alphabet_indexes[iterator] = 2;
				}
				else
				{
					if (((fm_index_entry[bitcounts[warp_index*32 + 3] + 9]) >> bitcounts[warp_index*32 + 2]) & 1)
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
				if (in_warp_index == 0)
				{
					bitcounts[tid] = indexes[warp_index*2 + 1] - bitcounts[warp_index*32];
				}
			}
			

			barrier(CLK_LOCAL_MEM_FENCE);


			//max 8 threadov, 9-17 pozicia
			//spocitat 'vlavo"
			if (in_warp_index >  8 && (alphabet_indexes[iterator] == 0 || alphabet_indexes[iterator] == 3))
			{

				if (in_warp_index<(9+bitcounts[warp_index*32]/32))
				{
					bitcounts[tid] = popcount(fm_index_entry[tid]);
				}
			
				else if (in_warp_index == (9+bitcounts[warp_index*32]/32))
				{	
					if (bitcounts[warp_index*32]%32){
						bitcounts[tid] = popcount(fm_index_entry[tid] >> (32 - bitcounts[warp_index*32]%32)); 
					}
					else
						bitcounts[tid] = 0;
				}
			}

			//spocitat "vpravo"
			if (in_warp_index >  8 && (alphabet_indexes[iterator] == 1 || alphabet_indexes[iterator] == 2))
			{

				a = fm_index_entry[warp_index*32 + 8];

				if (in_warp_index<(9+bitcounts[warp_index*32]/32))
				{
					bitcounts[tid] = popcount(fm_index_entry[tid + a]);
				}
			
				else if (in_warp_index == (9+bitcounts[warp_index*32]/32))
				{	
					if (bitcounts[warp_index*32]%32){
						bitcounts[tid] = popcount((fm_index_entry[tid + a]) >> (32 - bitcounts[warp_index*32]%32)); 
					}
					else 
						bitcounts[tid] = 0;
				}

				
			}

			barrier(CLK_LOCAL_MEM_FENCE);

			if (in_warp_index == 9)
			{
				if (alphabet_indexes[iterator] == 3 || alphabet_indexes[iterator] == 1)
				{
					for (a=1;a<=bitcounts[warp_index*32]/32;a++)
						bitcounts[tid] += bitcounts[tid + a];
				}	
				else
				{
					for (a=1;a<=bitcounts[warp_index*32]/32;a++)
						bitcounts[tid] += bitcounts[tid + a];
					bitcounts[tid] = bitcounts[warp_index*32] - bitcounts[tid];
				}
			}


			//add occ_counter in fm index entry
			if (in_warp_index == 0)
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

			if (in_warp_index == 9)
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
	



	if (in_warp_index ==0 ){
		//printf("vraciam %d %d\n",8*workgroupID, warp_index);

		output[8*workgroupID + warp_index] = 33;
	}
}
