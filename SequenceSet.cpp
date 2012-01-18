/*
 * SequenceSet.cpp
 */

#include <iostream>
#include <vector>
#include "Fasta.h"
#include "SequenceSet.h"
using namespace std;

///////////////////////////
// since data and onePointProbs are two-dimensional but implemented as a one-dimensional array, it is handy to have some access functions:
///////////////////////////
/**
 * reads an element from a two-dimensional array which is implemented
 * as one-dimensional
 */
template<class T> T readFromTwoDimArray(T* array, int pitch, int x, int y) {
	return array[y * pitch + x];
}
/**
 * writes an element into a two-dimensional array which is implemented
 * as one-dimensional
 */
template<class T> void writeIntoTwoDimArray(T* array, int pitch, int x, int y,
		T elem) {
	array[y * pitch + x] = elem;
}
/**
 * increase element in a two-dimensional array which is implemented
 * as one-dimensional
 */
template<class T> void incElementInTwoDimArray(T* array, int pitch, int x,
		int y) {
	array[y * pitch + x]++;
}
///////////////////////////


SequenceSet::SequenceSet(char* fastaFilePath) {
	const int numChars = NUMPROTEINCHARS;

	//read file content
	vector<vector<int> > fileContent = Fasta::readFastaFile(fastaFilePath);
	numSequences = fileContent.size();
	sequenceLength = fileContent.front().size();
	cout << "number of sequences: " << numSequences << endl << "sequence length: " << sequenceLength << endl << endl;

	//put data into data array
	data = new unsigned char[sequenceLength * numSequences];
	for (int i = 0; i < sequenceLength; i++)
		for (int sequence = 0; sequence < numSequences; sequence++)
			writeIntoTwoDimArray<unsigned char> (data, numSequences, sequence, i, fileContent.at(sequence).at(i));

	//compute one point probabilities
	onePointProbs = new float[sequenceLength * numChars];
	for (int i = 0; i < sequenceLength; i++) {
		for (int character = 0; character < numChars; character++)
			writeIntoTwoDimArray(onePointProbs, numChars, character, i, 0.0f);
		for (int seq = 0; seq < numSequences; seq++)
			incElementInTwoDimArray(onePointProbs, numChars, readFromTwoDimArray(data, numSequences, seq, i), i);
		for (int character = 0; character < numChars; character++)
			writeIntoTwoDimArray(onePointProbs, numChars, character, i, readFromTwoDimArray(onePointProbs, numChars, character, i)
					/ numSequences);
	}
}

SequenceSet::~SequenceSet() {
	delete[] data;
	delete[] onePointProbs;
}

int SequenceSet::getSequenceLength() {
	return sequenceLength;
}

int SequenceSet::getNumberOfSequences() {
	return numSequences;
}

float SequenceSet::getOnePointProb(int character, int i) {
	const int numChars = NUMPROTEINCHARS;
	return readFromTwoDimArray(onePointProbs, numChars, character, i);
}

float* SequenceSet::getOnePointProbs() {
	return onePointProbs;
}

unsigned char* SequenceSet::getData() {
	return data;
}

unsigned char SequenceSet::getData(int sequenceNo, int i) {
	return readFromTwoDimArray(data, numSequences, sequenceNo, i);
}
