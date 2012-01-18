////////////////////////////////////////////////////////////////
// GdI 3 - WS 11/12 - Praktikum 2
////////////////////////////////////////////////////////////////
// Bitte tragen Sie hier die Namen und Matrikelnummern Ihrer
// ein (maximal 3).
//
// Bsp.:
// Boris Baumstumpf 0999999133
// Axel Axt 12345678910
// Bruno Schneewittchen 66666336
////////////////////////////////////////////////////////////////

#include <iostream>
#include <CL/opencl.h>
#include <cstring>
#include <cmath>

#include <cstdlib>

#include <sys/time.h>

#include "MIcomputation.h"
#include "Fasta.h"
#include "oclutil.h"

#define LOOPCOUNT	1

void computeMI(bool cpu, bool oclgpu, SequenceSet& sequences, Matrix<float>& MI)
{
	//calculate MI
	cout << "computing MI" << endl;
	if (cpu)
		computeMIonCPU(sequences, MI);
	else//GPU
		computeMIonGPU(sequences, MI, oclgpu);
}

void computeMIonCPU(SequenceSet& sequences,	Matrix<float>& MI) {
	const int numChars = NUMPROTEINCHARS;
	const int sequenceLength = sequences.getSequenceLength();
	const int numSequences = sequences.getNumberOfSequences();

	const double epsilon=1e-6;
	
	timeval start, end;
	
	gettimeofday(&start, 0);
	for (int k = 0; k < LOOPCOUNT; k++) {
		//iterate over all column combinations
		for (int j = 0; j < sequenceLength; j++) {
			for (int i = 0; i <= j; i++) {
				//absolute number of occurrences of character pairs x,y: N_ij(x,y)
				int twoPointOccs[numChars][numChars];
				memset(twoPointOccs, 0, sizeof(twoPointOccs));
				//iterate through all sequences and compute two-point occurrences
				for (int seq = 0; seq < numSequences; seq++)
					twoPointOccs[sequences.getData(seq, i)][sequences.getData(seq, j)]++;
/*
				puts("===START===");
				for (int m=0; m<numChars; m++) {
					for (int n=0; n<numChars; n++)
						printf("%d %d: %d\n", m, n, twoPointOccs[m][n]);
					puts("");
				}
				puts("===STOP ===");
*/
				double MI_ij = 0;
				//sum over all x and y
				for (int x = 0; x < numChars; x++) {
					if (sequences.getOnePointProb(x, i) < epsilon)
						continue;
					for (int y = 0; y < numChars; y++) {
						if (sequences.getOnePointProb(y, j) < epsilon || twoPointOccs[x][y] == 0)
							continue;
						double p_ij_xy = double(twoPointOccs[x][y]) / double(numSequences);
						MI_ij += p_ij_xy * log2(p_ij_xy / (sequences.getOnePointProb(x, i) * sequences.getOnePointProb(y, j)));
					}
				}
				MI.set(i, j, MI_ij);
			}
		}
	}
	gettimeofday(&end, 0);
	std::cout << "execution time: " 
			<< (end.tv_sec - start.tv_sec ) * 1000 +  ( end.tv_usec - start.tv_usec) / 1000
			<< " milliseconds" << std::endl;
}

void computeMIonGPU(SequenceSet& sequence, Matrix<float>& MI, bool GPU)
{
	// initializes context and kernel and stores them
	OCL ocl(GPU);

	cl_int oclError1, oclError2;

	timeval start, end;

	// memory sizes
	size_t sequenceLength = sequence.getSequenceLength();
	size_t numSequences = sequence.getNumberOfSequences();
	
	// matrix MI is of size numElements
	size_t numElements = sequenceLength * sequenceLength;
	size_t sequenceSize = sequence.getNumberOfSequences() * sequenceLength;
	size_t onePointProbsSize = sequenceLength * NUMPROTEINCHARS;
		
	// host memory
	float * dst = new float[MI.size()];
	memset(dst, 0, MI.size());
	
	// device memory for sequences, one point probablities and resulting matrix
	cl_mem oclDevSrcSequence, oclDevSrcOnePointProbs, oclDevDstMI;
	
	// size for a work group: each workgroup computes one matrix entry, thus computes the correlation
	// one time for each character => 25 work items are sufficient
	size_t localWorkSize[2] = { 1, 1 };
	
	// global work size defines the total amount of threads over all work group, thus needs to be a multiple of the local
	// work size in each dimension.
	size_t globalWorkSize[2] = { localWorkSize[0] * sequenceLength, localWorkSize[1] * sequenceLength };
	
	// create buffer on device, one for each input array
	oclDevSrcSequence = clCreateBuffer(		ocl.oclContext,
											CL_MEM_READ_ONLY,
											sizeof(cl_uchar) * sequenceSize,
											0, &oclError1);

	oclDevSrcOnePointProbs = clCreateBuffer(ocl.oclContext,
											CL_MEM_READ_ONLY,
											sizeof(cl_float) * onePointProbsSize,
											0, &oclError2);
	oclError1 |= oclError2;

	oclDevDstMI = clCreateBuffer(			ocl.oclContext,
											CL_MEM_WRITE_ONLY,
											sizeof(cl_float) * numElements,
											0, &oclError2);
	oclError1 |= oclError2;
	
	if (oclError1 != CL_SUCCESS) {
		std::cout << "error while allocating buffers" << std::endl;
		exit(1);
	}
	
	// set buffer to appropriate kernel arguments
	oclError1 = clSetKernelArg(ocl.oclKernel, 0, sizeof(cl_mem), (void*)&oclDevSrcSequence);
	oclError1 |= clSetKernelArg(ocl.oclKernel, 1, sizeof(cl_mem), (void*)&oclDevSrcOnePointProbs);
	oclError1 |= clSetKernelArg(ocl.oclKernel, 2, sizeof(cl_mem), (void*)&oclDevDstMI);
	oclError1 |= clSetKernelArg(ocl.oclKernel, 3, sizeof(cl_uint), &sequenceLength);
	oclError1 |= clSetKernelArg(ocl.oclKernel, 4, sizeof(cl_uint), &numSequences);

	if (oclError1 != CL_SUCCESS) {
		std::cout << "error while setting arguments: " << ocl.oclErrorString(oclError1) << std::endl;
		exit(1);
	}

	// copy host memory to device, non-blocking copy
	oclError1 = clEnqueueWriteBuffer(	ocl.oclCmdQueue,
						oclDevSrcSequence,
						CL_FALSE,
						0,
						sizeof(cl_uchar) * sequenceSize,
						(const void *) sequence.getData(),
						0, 0, 0);

	oclError1 |= clEnqueueWriteBuffer(	ocl.oclCmdQueue,
						oclDevSrcOnePointProbs,
						CL_FALSE,
						0,
						sizeof(cl_float) * onePointProbsSize,
						(const void *) sequence.getOnePointProbs(),
						0, 0, 0);

	if (oclError1 != CL_SUCCESS) {
		std::cout << "error while writing to device " << ocl.oclErrorString(oclError1) << std::endl;
		exit(1);
	}

	// execute kernel LOOPCOUNT times and measure execution time
	// TODO LOOPCOUNT aendern, um Kernel mehrfach auszufuehren
	gettimeofday(&start, 0);
	for (int i = 0; i < LOOPCOUNT; ++i) {
		oclError1 = clEnqueueNDRangeKernel(	ocl.oclCmdQueue,
							ocl.oclKernel,
							2,			// dimension
							0,
							globalWorkSize,
							localWorkSize,
							0, 0, 0);
	
	
		if (oclError1 != CL_SUCCESS) {
			std::cout << "error while executing kernel: " << ocl.oclErrorString(oclError1) << std::endl;
			exit(1);
		}
	}
	
	// clFinish blocks until all issued commands so far are completed, necessary for computing execution time
	oclError1 = clFinish(ocl.oclCmdQueue);
	gettimeofday(&end, 0);

	// read memory from device, store in temporary array and if no error happend copy to result matrix
	oclError1 = clEnqueueReadBuffer(	ocl.oclCmdQueue,
						oclDevDstMI,
						CL_TRUE,
						0,
						sizeof(cl_float) * numElements,
						dst,
						0, 0, 0);

	if (oclError1 != CL_SUCCESS) {
		std::cout << "error while reading from device: " << ocl.oclErrorString(oclError1) << std::endl;
		exit(1);
	}

	std::cout << "execution time: " 
			<< (end.tv_sec - start.tv_sec ) * 1000 +  ( end.tv_usec - start.tv_usec) / 1000
			<< " milliseconds" << std::endl;
	
	// fill the matrix with the computed results
	MI.copyElements(dst);
	
	// release used memory, can cause really bad crashes otherwise
	clReleaseMemObject(oclDevSrcSequence);
	clReleaseMemObject(oclDevSrcOnePointProbs);
	clReleaseMemObject(oclDevDstMI);
}
