/*
 * MIcomputation.cuh
 *
 * Performs the whole MI computation. Dispatches the task either to the CPU or the GPU depending on the given parameters.
 *
 */

#ifndef MICOMPUTATION_H_
#define MICOMPUTATION_H_

#include "Matrix.h"
#include "SequenceSet.h"

/*
 * performs the MI computation
 * cpu: whether the computation shall be performed on the cpu or on the gpu
 * sequences: the sequence set to operate on
 * MI: the resulting MI will be stored in this matrix (be careful to define MI with the right dimensions)
 */
void computeMI(bool cpu, bool oclgpu, SequenceSet& sequences, Matrix<float>& MI);
/*
 * Calculates MI on the CPU.
 */
void computeMIonCPU(SequenceSet& sequences, Matrix<float>& MI);
/*
 * Calculates MI on the GPU.
 */
void computeMIonGPU(SequenceSet& sequences, Matrix<float>& MI, bool GPU);

#endif /* MICOMPUTATION_H_ */
