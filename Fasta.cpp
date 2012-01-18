/*
 * Fasta.cpp
 */

#include "Fasta.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>

vector<vector<int> > Fasta::readFastaFile(char* filePath) {
	vector<vector<int> > result;
	string invalidChars;
	string line;
	ifstream myfile(filePath);
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			getline(myfile, line);
			if (line.length() == 0)
				continue;
			//if '>' occurs, a new protein sequence starts in the next line
			if (line.at(0) == '>') {
                                result.push_back(vector<int>());
				continue;
			}
			//remove windows style carriage return
			size_t carrRetPos = line.find_first_of('\r');
			if (carrRetPos != string::npos)
                                line.replace(carrRetPos, 1, "");

                        //convert data and put it into data array (since no new '>' has been
                        //found, no new sequence has started and no pushback is nescessary)
                        for (unsigned int i = 0; i < line.length(); i++)
                                result.back().push_back(Fasta::charToNumber(line.at(i), invalidChars));

		}
		myfile.close();

		//check whether all sequences are of equal length
                for (unsigned int i = 1; i < result.size(); i++)
                        if (result.at(i).size() != result.front().size()) {
				cerr << "length of sequence no. " << i + 1 << " doesn't match with the length of the 1st sequence." << endl;
				exit(-1);
			}
	} else {
		cerr << "Unable to open file" << endl;
		exit(-1);
	}


	if (!invalidChars.empty())
		cout << "The file contained invalid characters: " << invalidChars << "\nAll invalid characters have been changed to \'X\'." << endl << endl;
	return result;
}

const string Fasta::proteinChars = "ACDEFGHIKLMNPQRSTVWYX-";

char Fasta::charToNumber(char c, string& invalidChars) {
	//set to upper case
	if (c >= 'a' && c <= 'z')
		c -= ('a' - 'A');

	string validChars;
	validChars = proteinChars;
	
	size_t pos = validChars.find(c);
	if (pos == string::npos) {
		if (invalidChars.find(c) == string::npos)
			invalidChars.push_back(c);
		return Fasta::charToNumber('X', invalidChars);
	} else
		return pos;
}
