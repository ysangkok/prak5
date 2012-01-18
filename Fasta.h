/*
 * Fasta.h
 *
 * This class serves to read a fasta file and
 * convert the fasta characters into numbers.
 *
 */

#include <vector>
#include <string>

#ifndef FASTA_H_
#define FASTA_H_

#define NUMPROTEINCHARS 22 //"ACDEFGHIKLMNPQRSTVWYX-" consists of 22 characters
#define PROTEINNUMCHARSHALFNEXTPOWEROFTWO 16//16 is the next power of two from 22/2
using namespace std;

class Fasta {
public:
	/**
	 * reads a fasta file
	 * every string in the resulting vector of strings is
	 * equivalent to one protein/dna sequence from the file
	 */
	static vector<vector<int> > readFastaFile(char* filePath);

private:
	/**
	 * translates a character from a FASTA file into the
	 * equivalent number (e.g. 'A' --> 0, '-' --> 21)
	 *
	 * if the character is illegal, it will be put into the
	 * invalidChars string (if it's not already in there)
	 */
	static char charToNumber(char c, string& invalidChars);
	static const string proteinChars;
};

#endif /* FASTA_H_ */
