#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "KmerChaosGen.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2019 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

const char progName[] = "buildkmerChaos";

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName, int& kmerSize);

int main(int argc,char *argv[]){

	string inFileName;
	string outFileName;
	int kmerSize = -1;
	
	if(!getInputs(argc, argv, inFileName, outFileName, kmerSize)){
		cerr << "Process aborted.\n";
		return 0;
	}
	
	KmerChaosGen kmerChaosGen(inFileName, outFileName, kmerSize);
	
	if(kmerChaosGen.GenKmerChaos()){
		return 0;
	}else{
		return 1;
	}
}

bool getInputs(int argc, char* argv[], string& inFileName, string& outFileName, int& kmerSize){
	if(argc != 4){
		cerr << "\t***** " << progName << " *****\n\t- Andrew Spriggs, CSIRO Ag&Food, 2019 -\n";
		cerr << "Build chaos game representation of kmers of sequence, as per:\n";
		cerr << "https://towardsdatascience.com/chaos-game-representation-of-a-genetic-sequence-4681f1a67e14.\n";
		cerr << "Input file may be fasta or fastq and may be .gz compressed.\n";
		cerr << "Correct command line usage:\n" << argv[0] << " <infile> <outfile> <kmer size>\n";
		return false;
	}
	
	inFileName = argv[1];
	outFileName = argv[2];
	kmerSize = atoi(argv[3]);
	
	return true;
}

