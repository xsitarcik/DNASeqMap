#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <CL/cl.h>


#define MEM_SIZE 128
#define MAX_SOURCE_SIZE 0x100000
#define DEVICE_GROUP_SIZE 256

void resize_array(int** orig, int newSize) {
	printf("Resizing to: %lu\n", sizeof(*orig) * newSize);
	printf("before size: %lu\n",sizeof(*orig));
	int *temp = realloc(*orig, sizeof(*orig) * newSize);

	if(temp == NULL) {
		puts("resize failed");
		exit(EXIT_FAILURE);
	} else {
		*orig = temp;
	}
}

void error_handler(char err[], int code) {
	if(code != CL_SUCCESS) {
		printf("%s, Error code:%d\n", err, code);
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char * argv[]) {
	int ret;
	int orig_inputSize = atoi(argv[1]);
	int inputSize,residue;
	printf("Specified input array will have %d items\n", orig_inputSize);
	printf("The program will sum the integers in the range [0, %d)\n", orig_inputSize);
	printf("==============\n");

	/* create input data */
	int* input = (int *)malloc(sizeof(*input) * orig_inputSize);
	for(int i = 0; i < orig_inputSize; i++) {
		input[i] = 1;
	}

	printf("device group size je %d\n",DEVICE_GROUP_SIZE);
	int sizeDiff = orig_inputSize - DEVICE_GROUP_SIZE;
	int multiplier = 1;
	
	if(sizeDiff < 0) {
		resize_array(&input, DEVICE_GROUP_SIZE);
	}
	else if(sizeDiff > 0) {
		multiplier = ((orig_inputSize % (DEVICE_GROUP_SIZE*2)) == 0) ? (orig_inputSize / (DEVICE_GROUP_SIZE*2)): (orig_inputSize / (DEVICE_GROUP_SIZE*2)) + 1;
			/* read results from output buffer */
		multiplier = multiplier*2;
		resize_array(&input, multiplier * DEVICE_GROUP_SIZE);
		//multiplier = multiplier/2;
	}

	for(int i = 0; i < abs(sizeDiff); i++) {
		input[orig_inputSize + i] = 0;
	}

	

	printf("mulitplier %d * device group %d = %d\n",multiplier,DEVICE_GROUP_SIZE, multiplier*DEVICE_GROUP_SIZE);
	fflush(stdout);
	inputSize = multiplier * DEVICE_GROUP_SIZE;

	printf("multiplier je %d\n",multiplier);

	if (multiplier%2)
		residue = DEVICE_GROUP_SIZE - (inputSize-orig_inputSize);
	else
		residue = 0;
	printf("residue je %d\n",residue);

	multiplier = (multiplier)/2;


	printf("NICE\n");
	fflush(stdout);
	/* define platform */
	cl_platform_id platformID;
	ret = clGetPlatformIDs(1, &platformID, NULL);
	error_handler("Get platform error", ret);	

	printf("NICE\n");
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

	printf("NICE\n");
	fflush(stdout);
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
	inputMemObj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int)*inputSize, NULL, &ret);
	error_handler("Create input buffer failed", ret);
	ret = clEnqueueWriteBuffer(commandQueue, inputMemObj, CL_TRUE, 0, sizeof(int)*inputSize, (const void*)input, 0, NULL, NULL);
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
	
	printf("MERDE\n");
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
	ret = clSetKernelArg(kernel, 2, sizeof(int), (void *)&inputSize);
	error_handler("Set arg 3 failure", ret);
	ret = clSetKernelArg(kernel, 3, sizeof(int) * DEVICE_GROUP_SIZE, NULL);
	error_handler("Set arg 4 failure", ret);

	/* enqueue and execute */
	const size_t globalWorkSize = inputSize/2;
	const size_t localWorkSize = DEVICE_GROUP_SIZE;
	printf("global %d, local %d\n",globalWorkSize, localWorkSize);
	ret = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL);
	error_handler("Enqueue/execute failure", ret);
	

	clock_t begin = clock();

	int * results = (int*) malloc(sizeof(int)*multiplier);
	ret = clEnqueueReadBuffer(commandQueue, outputMemObj, CL_TRUE, 0, sizeof(int)*multiplier, results, 0, NULL, NULL);
	error_handler("Read output buffer fail", ret);
	ret = clEnqueueReadBuffer(commandQueue, outputMemObj, CL_TRUE, 0, sizeof(int)*multiplier, results, 0, NULL, NULL);
	error_handler("Read output buffer fail", ret);

	//je nutne spocitat residue
	int sum = 0;
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


	/* release */
	ret = clReleaseKernel(kernel);
	ret = clReleaseProgram(program);
	ret = clReleaseMemObject(inputMemObj);
	ret = clReleaseMemObject(outputMemObj);
	ret = clReleaseCommandQueue(commandQueue);
	ret = clReleaseContext(context);	

	free(input);
	free(results);

	exit(EXIT_SUCCESS);
}