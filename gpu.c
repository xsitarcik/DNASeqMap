#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <CL/cl.h>


#define MEM_SIZE 128
#define MAX_SOURCE_SIZE 0x100000
#define DEVICE_GROUP_SIZE 256
#define NUM_OF_READS_PER_WARP 4 

/*void resize_array(int** orig, int newSize) {
	printf("Resizing to: %lu\n", sizeof(*orig) * newSize);
	printf("before size: %lu\n",sizeof(*orig));
	int *temp = realloc(*orig, sizeof(*orig) * newSize);

	if(temp == NULL) {
		puts("resize failed");
		exit(EXIT_FAILURE);
	} else {
		*orig = temp;
	}
}*/

void error_handler(char err[], int code) {
	if(code != CL_SUCCESS) {
		printf("%s, Error code:%d\n", err, code);
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char * argv[]) {
	int ret;
 	int ALPHABET_SIZE = 4;
 	int PATTERN_LENGTH = 64;

	unsigned int * count_table = (unsigned int *)malloc(sizeof(unsigned int) *  (ALPHABET_SIZE+1));
	count_table[0] = 4;
	count_table[1] = 5;
	count_table[2] = 6;
	count_table[3] = 7;
	count_table[4] = 9;

	/* create input data */
	int orig_inputSize = 1024;
	int multiplier = 4;

	char* input = (char *)malloc(sizeof(char) *  1500);
	input = "TCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCATCTATGGAAACTACAGGACTAACCTTCCTGGCAACCGGGGGCTGGGAATCTGTCACATGAGTCA";

	int inputSize = 1024;

	


	printf("%d %d %d %d\n",count_table[0],count_table[1],count_table[2],count_table[3]);
	fflush(stdout);
	/* define device */
	cl_device_id deviceID;
	ret = clGetDeviceIDs(platformID, CL_DEVICE_TYPE_GPU, 1, &deviceID, NULL);
	error_handler("Get deviceID error", ret);

	cl_char vendorName[1024] = {0};
	cl_char deviceName[1024] = {0};
	ret = clGetDeviceInfo(deviceID, CL_DEVICE_VENDOR, sizeof(vendorName), vendorName, NULL);
	ret = clGetDeviceInfo(deviceID, CL_DEVICE_NAME, sizeof(deviceName), deviceName, NULL);
	printf("Connecting to %s %s...\n", vendorName, deviceName);

	/* define context */
	cl_context context;
	context = clCreateContext(NULL, 1, &deviceID, NULL, NULL, &ret);
	error_handler("define context error", ret);

	/* define command queue */
	cl_command_queue commandQueue;
	commandQueue = clCreateCommandQueue(context, deviceID, 0, &ret);
	error_handler("Create command queue error", ret);

	/* create memory objects */
	cl_mem inputMemObj;
	inputMemObj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(char)*inputSize, NULL, &ret);
	error_handler("Create input buffer failed", ret);
	ret = clEnqueueWriteBuffer(commandQueue, inputMemObj, CL_TRUE, 0, sizeof(char)*inputSize, (const void*)input, 0, NULL, NULL);
	error_handler("Write to input buffer failed", ret);

	/* create memory objects */
	cl_mem inputMemObj_count_table;
	inputMemObj_count_table = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned int)*(ALPHABET_SIZE+1), NULL, &ret);
	error_handler("Create input buffer failed", ret);
	ret = clEnqueueWriteBuffer(commandQueue, inputMemObj_count_table, CL_TRUE, 0, sizeof(unsigned int)*(ALPHABET_SIZE+1), (const void*)count_table, 0, NULL, NULL);
	error_handler("Write to input buffer failed", ret);

	cl_mem inputMemObj_fm_index;
	inputMemObj_fm_index = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned int)*(ALPHABET_SIZE+1), NULL, &ret);
	error_handler("Create input buffer failed", ret);
	ret = clEnqueueWriteBuffer(commandQueue, inputMemObj_fm_index, CL_TRUE, 0, sizeof(unsigned int)*(ALPHABET_SIZE+1), (void*)count_table, 0, NULL, NULL);
	error_handler("Write to input buffer failed", ret);

	// output need only be an array of inputSize / DEVICE_GROUP_SIZE long
	cl_mem outputMemObj;
	

	outputMemObj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int)*multiplier, NULL, &ret);
	error_handler("Create output buffer failed", ret);

	/* read kernel file from source */
	char fileName[] = "./sum.cl";
	char *sourceStr;	
	size_t sourceSize;

	FILE *fp = fopen(fileName, "r");
	if(!fp) {
		puts("Failed to load kernel file");
		exit(EXIT_FAILURE);
	}
	sourceStr = (char*)malloc(MAX_SOURCE_SIZE);
	sourceSize = fread(sourceStr, 1, MAX_SOURCE_SIZE, fp);
	fclose(fp);


	/* create program object */
	cl_program program = clCreateProgramWithSource(context, 1, (const char**)&sourceStr, (const size_t*)&sourceSize, &ret);	
	error_handler("create program failure", ret);
	ret = clBuildProgram(program, 1, &deviceID, NULL, NULL, NULL);
	if(ret != CL_SUCCESS) {
		puts("Build program error");
		size_t len;
		char buildLog[2048];
		clGetProgramBuildInfo(program, deviceID, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, &len);
		printf("%s\n", buildLog);
	}

	/* create kernel */
	cl_kernel kernel = clCreateKernel(program, "sum", &ret);
	error_handler("Create kernel failure", ret);
	
	/* set kernel arguments */
	ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&inputMemObj);
	error_handler("Set arg 1 failure", ret);
	ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&outputMemObj);
	error_handler("Set arg 2 failure", ret);
	ret = clSetKernelArg(kernel, 2, sizeof(int), (void *)&PATTERN_LENGTH);
	error_handler("Set arg 3 failure", ret);
	ret = clSetKernelArg(kernel, 3, sizeof(unsigned int) * 2 * NUM_OF_READS_PER_WARP, NULL);
	error_handler("Set arg 4 failure", ret);
	ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&inputMemObj_count_table);
	error_handler("Set arg 5 failure", ret);
	ret = clSetKernelArg(kernel, 5, sizeof(unsigned int) * 5, NULL);
	error_handler("Set arg 6 failure", ret);
	ret = clSetKernelArg(kernel, 6, sizeof(char) * NUM_OF_READS_PER_WARP * PATTERN_LENGTH, NULL);
	error_handler("Set arg 7 failure", ret);
	ret = clSetKernelArg(kernel, 7, sizeof(unsigned int) * 4 * NUM_OF_READS_PER_WARP, NULL); //indexes
	error_handler("Set arg 9 failure", ret);
	ret = clSetKernelArg(kernel, 8, sizeof(unsigned int) * DEVICE_GROUP_SIZE, NULL); //fm_index entry
	error_handler("Set arg 10 failure", ret);
	ret = clSetKernelArg(kernel, 9, sizeof(unsigned int) * DEVICE_GROUP_SIZE * 2, NULL); //count_table_results
	error_handler("Set arg 11 failure", ret);
	ret = clSetKernelArg(kernel, 10, sizeof(unsigned int) * DEVICE_GROUP_SIZE, NULL); //bitcounts
	error_handler("Set arg 12 failure", ret);
	ret = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *)&inputMemObj_fm_index);
	error_handler("Set arg 5 failure", ret);


	/* enqueue and execute */
	const size_t globalWorkSize = 1024/2;
	const size_t localWorkSize = 256;
	printf("global %d, local %d\n",globalWorkSize, localWorkSize);
	ret = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL);
	error_handler("Enqueue/execute failure", ret);
	

	//clock_t begin = clock();

	int * results = (int*) malloc(sizeof(int)*multiplier);

	ret = clEnqueueReadBuffer(commandQueue, outputMemObj, CL_TRUE, 0, sizeof(int)*multiplier, results, 0, NULL, NULL);
	error_handler("Read output buffer fail", ret);

	int sum = 0;
	//printf("multiplier je %d\n",multiplier);
	for(int i = 0; i < multiplier; i++) {
		printf("wat %d\n",results[i]);
		sum += results[i];
	}

	/*ret = clEnqueueReadBuffer(commandQueue, outputMemObj, CL_TRUE, 0, sizeof(int)*multiplier, results, 0, NULL, NULL);
	error_handler("Read output buffer fail", ret);*/

	//je nutne spocitat residue
	/*int sum = 0;
	for (int i = orig_inputSize-residue;i<orig_inputSize;i++){
		printf("%d\n", input[i]);
		sum += input[i]; 
	}

	
	//printf("multiplier je %d\n",multiplier);
	for(int i = 0; i < multiplier; i++) {
		//printf("wat %d\n",results[i]);
		sum += results[i];
	}
	//printf("posledne je %d\n",results[multiplier]);
	clock_t end = clock();


 double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
 printf("construction of wavelet tree took %lf seconds\n", time_spent);
 fflush(stdout);

	printf("==============\n");
	printf("The sum of the array is: %d\n", sum);
		
		clock_t begin2 = clock();
		sum = 0;
	for(int i = 0; i < inputSize; i++) {
		sum += input[i];
	}
	clock_t end2 = clock();
 double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
 printf("construction of wavelet tree took %lf seconds\n", time_spent2);
 printf("The sum of the array is: %d\n", sum);
 fflush(stdout);
*/

	/* release */
	ret = clReleaseKernel(kernel);
	ret = clReleaseProgram(program);
	ret = clReleaseMemObject(inputMemObj);
	ret = clReleaseMemObject(outputMemObj);
	ret = clReleaseCommandQueue(commandQueue);
	ret = clReleaseContext(context);	

	printf("MEKO\n");
	fflush(stdout);
	//free(input);
	free(results);

	exit(EXIT_SUCCESS);
}