__kernel void sum(__global const int *input, __global int *output, const int inputSize, __local int* reductionSums) {
	const int tid = get_local_id(0); //tid
	//const int globalID = BLOCK_SIZE; //i
	const int blockDim = get_local_size(0);//blockDim.x
	const int workgroupID = get_group_id(0); //blockIdx.x
	int i = workgroupID*blockDim*2+tid;
	const int gridSize = get_num_groups(0)*2*blockDim;
	//const int i = workgroupID*blockDim+tid;
	reductionSums[tid] = input[i]+input[i+blockDim]; //sdata[tid]=g_data[i]
	/*while (i<inputSize){
		reductionSums[tid] = input[i]+input[i+blockDim]; //sdata[tid]=g_data[i]

		i += gridSize;
	}*/
	
	//reductionSums[tid] = input[i];

	/*for(int offset = blockDim / 2; offset > 0; offset >>= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);	// wait for all other work-items to finish previous iteration.
		
		if(tid < offset) {
			reductionSums[tid] += reductionSums[tid + offset];

		}
	}*/

	if (tid<128){
		barrier(CLK_LOCAL_MEM_FENCE);	// wait for all other work-items to finish previous iteration.
		reductionSums[tid] += reductionSums[tid + 128];
	}	

	if (tid<64){
		barrier(CLK_LOCAL_MEM_FENCE);	// wait for all other work-items to finish previous iteration.
		reductionSums[tid] += reductionSums[tid + 64];
	}
	if (tid<32){
		barrier(CLK_LOCAL_MEM_FENCE);	// wait for all other work-items to finish previous iteration.
		reductionSums[tid] += reductionSums[tid + 32];
		barrier(CLK_LOCAL_MEM_FENCE);	// wait for all other work-items to finish previous iteration.
		reductionSums[tid] += reductionSums[tid + 16];
		barrier(CLK_LOCAL_MEM_FENCE);	// wait for all other work-items to finish previous iteration.
		reductionSums[tid] += reductionSums[tid + 8];
		barrier(CLK_LOCAL_MEM_FENCE);	// wait for all other work-items to finish previous iteration.
		reductionSums[tid] += reductionSums[tid + 4];
		barrier(CLK_LOCAL_MEM_FENCE);	// wait for all other work-items to finish previous iteration.
		reductionSums[tid] += reductionSums[tid + 2];
		barrier(CLK_LOCAL_MEM_FENCE);	// wait for all other work-items to finish previous iteration.
		reductionSums[tid] += reductionSums[tid + 1];

	}	


	/*if (tid < 32){
		//barrier(CLK_LOCAL_MEM_FENCE);
		reductionSums[tid] += reductionSums[tid + 32];
		//barrier(CLK_LOCAL_MEM_FENCE);
		reductionSums[tid] += reductionSums[tid + 16];
		//barrier(CLK_LOCAL_MEM_FENCE);
		reductionSums[tid] += reductionSums[tid + 8];
		//barrier(CLK_LOCAL_MEM_FENCE);
		reductionSums[tid] += reductionSums[tid + 4];
		//barrier(CLK_LOCAL_MEM_FENCE);
		reductionSums[tid] += reductionSums[tid + 2];
		//barrier(CLK_LOCAL_MEM_FENCE);
		reductionSums[tid] += reductionSums[tid + 1];
	}*/
	if(tid == 0) {	// the root of the reduction subtree
		output[workgroupID] = reductionSums[0];
	}
}

__kernel void inefficient(__global const int *input, __global int *output, const int inputSize, __local int* reductionSums) {
	const int globalID = get_global_id(0);
	const int localID = get_local_id(0);
	const int localSize = get_local_size(0);
	const int workgroupID = globalID / localSize;

	reductionSums[localID] = input[globalID];
	
	barrier(CLK_LOCAL_MEM_FENCE); // wait for the rest of the work items to copy the input value to their local memory.
	if(localID == 0) {
		int sum = 0;
		for(int i = 0; i < localSize; i++) {
			if((i+localID) < localSize) {
				sum += reductionSums[i+localID];
			}
		}
		output[workgroupID] = sum;
	}
}

