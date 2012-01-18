/*
 * SequenceSet.h
 *
 */

#ifndef SEQUENCESET_H_
#define SEQUENCESET_H_

class SequenceSet {
public:
	//Reads a set of sequences from a FASTA file.
	SequenceSet(char* fastaFilePath);
	virtual ~SequenceSet();
	int getSequenceLength();
	int getNumberOfSequences();
	float getOnePointProb(int character, int i);
	float* getOnePointProbs();
	unsigned char* getData();
	unsigned char getData(int sequenceNo, int i);
private:
	//number of sequences
	int numSequences;
	//length of the sequences
	int sequenceLength;
	/**
	 * the sequences
	 * NOTE: this is a two-dimensional array in a one-dimensional
	 */
	unsigned char* data;
	/**
	 * the one-point probabilities p_i(x) that a character x appears in
	 * position i of any of the sequences
	 * NOTE: this is a two-dimensional array in a one-dimensional
	 */
	float* onePointProbs;
};

#endif /* SEQUENCESET_H_ */
