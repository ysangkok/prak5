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

#define NUMCHARS 	22		// amount of possible characters in a sequence
#define SEQLENGTH 	1024	// maximum sequence length
#define NUMSEQ		512		// maximum number of sequences

//#pragma OPENCL EXTENSION cl_amd_printf : enable

__kernel void oclMIcomputation(__global uchar * sequences, __global float * onePointProbs, __global float * result, uint sequenceLength, uint numSequences)
{
	const float epsilon = 1e-6f;

	const int x = get_global_id(0), y = get_global_id(1), xWid = get_global_size(0);
	if (!(x <= y)) return;
	//printf("%d, %d, %d, %d\n", get_global_id(0), get_global_id(1), get_local_id(0), get_local_id(1));


	int twoPointOccs[NUMCHARS][NUMCHARS];
	int i;
	for (i = 0; i<NUMCHARS*NUMCHARS; i++) ((int*) twoPointOccs)[i]=0;

	for (int seq = 0; seq < numSequences; seq++)
		twoPointOccs[sequences[x * numSequences + seq]][sequences[y * numSequences + seq]]++;

/*
	printf("===START===\n");
	for (int m=0; m<NUMCHARS; m++) {
		for (int n=0; n<NUMCHARS; n++)
			printf("%d %d: %d\n", m, n, twoPointOccs[m][n]);
		printf("\n");
	}
	printf("===STOP ===\n");
*/

	float MI_ij = 0;
	for (int x1 = 0; x1 < NUMCHARS; x1++) {
		if (onePointProbs[x*NUMCHARS+x1] < epsilon)
			continue;
		for (int y1 = 0; y1 < NUMCHARS; y1++) {
			if (onePointProbs[y*NUMCHARS+y1] < epsilon || twoPointOccs[x1][y1] == 0)
				continue;
			float p_ij_xy = ((float) (twoPointOccs[x1][y1])) / ((float) numSequences);
			MI_ij += p_ij_xy * log2(p_ij_xy / (onePointProbs[x*NUMCHARS+x1] * onePointProbs[y*NUMCHARS+y1]));
		}
	}

	result[y*xWid + x] = MI_ij;
}
