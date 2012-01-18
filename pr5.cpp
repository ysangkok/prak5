#include <iostream>
#include <fstream>
#include <vector>
#include "Fasta.h"
#include "SequenceSet.h"
#include "Matrix.h"
#include "MIcomputation.h"
using namespace std;

int main(int argc, char *argv[]) {

	//parameter parsing
	vector < string > params;
	for (int i = 0; i < argc; i++)
		params.push_back(string(argv[i]));
	string inpath, MI_outpath;
	bool cpu = false, oclgpu = false, fast = false;
	const string instring = "-i=", cpustring = "-cpu", uppertrianglestring = "-upperTriangle", oclgpustring = "-oclgpu", fststr="-fast";
	for (unsigned i = 0; i < params.size(); i++)
		if (params[i].find(instring) != string::npos)
			inpath = params[i].erase(0, instring.length());
		else if (params[i].find(cpustring) != string::npos)
			cpu = true;
		else if (params[i].find(oclgpustring) != string::npos)
			oclgpu = true;
		else if (params[i].find(fststr) != string::npos)
			fast = true;

	if (inpath.empty()) {
		cerr << "Error: no input file has been provided\n" << "parameter syntax is:\n"
				<< instring << "<infile>\n[" 
				<< cpustring << " for CPU computation (default: GPU)]\n["
				<< oclgpustring << " use GPU to compute via OpenCL, is ignored if " << cpustring << " is on (default: OCLCPU)]"
				<< endl;
		return -1;
	}
	MI_outpath = string(inpath).append(".MI.out");
	cout << "inpath: " << inpath << endl << "MI outpath: " << MI_outpath << endl;
	cout << "computing on the " << (cpu ? "CPU" : (oclgpu ? "GPU" : "CPU with OpenCL")) << endl << endl;
	remove(MI_outpath.data());

	//sequence reading
	SequenceSet sequences((char*) inpath.data());

	//shuffle null model computation
	int len = sequences.getSequenceLength();
	Matrix<float> MI(len, len);
	//time_t timer = time(NULL);
	//cout << "starting computation at: " << ctime(&timer) << endl;
	computeMI(cpu, oclgpu, sequences, MI);
	//timer = time(NULL);
	//cout << endl << "computation finished at: " << ctime(&timer) << endl;

	//output
	ofstream myfile;
	cout << "writing output to " << MI_outpath << endl;
	myfile.open((char*) MI_outpath.data(), ios_base::app);
	if (fast) {
		MI.fastprint(myfile, true);
	} else {
		MI.print(myfile, true);
	}
	myfile.close();
}
