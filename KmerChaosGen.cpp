#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctype.h>
#include <math.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "KmerChaosGen.h"
#include "SeqReader.h"
using namespace boost::iostreams;

/* By Andrew Spriggs, CSIRO Ag&Food, 2019 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

/*** Initialise
**/
KmerChaosGen::KmerChaosGen(const string& aInFileName, const string& aOutFileName, const int akmerSize ){
	prepareKmerChaosGen(aInFileName, aOutFileName, akmerSize);
	return;
}

/*** Actual constructor 
**/
void KmerChaosGen::prepareKmerChaosGen(const string& aInFileName, const string& aOutFileName, const int akmerSize){

	ready = false;
	inFileName = aInFileName;
	outFileName = aOutFileName;
	kmerSize = akmerSize;

	// Check kmer size requested
	if(kmerSize < 2 || kmerSize > 8){
		cerr << "Invalid value for kmer size!\nNo processing will follow.\n";
		return;
	}
	
	// Check that output file is openable
	outFile.open(aOutFileName.c_str());
	if(!outFile.is_open()){
		cerr << "Unable to open output file " << aOutFileName << "!\nNo processing will follow.\n";
		return;
	}

	// Initialise matrix to zeroes
	quadrantWidth = pow(2.0, kmerSize)/2;
	for(int i=0; i<quadrantWidth; i++){
		MatrixRow aRow(quadrantWidth, 0);
		quadrantA.push_back(aRow);
		quadrantT.push_back(aRow);
		quadrantC.push_back(aRow);
		quadrantG.push_back(aRow);
	}

	ready = true;
	return;
}


KmerChaosGen::~KmerChaosGen(){
	if(outFile.is_open()){
		outFile.close();
	}
}

/*** Launch the full kmer analysis process
**/
bool KmerChaosGen::GenKmerChaos(){
	
	if(!ready){
		return false;
	}

	SeqReader inFile(inFileName);
	while(inFile.nextSeq()){
		int sLen = inFile.getSeqLen();
		if(sLen >= kmerSize){
			string seq = inFile.getSeq();
			#pragma omp parallel for
			for(int i=0; i < 4; i++){
				switch(i){
					case 0:
						processSeqInQuadrant(seq, sLen, quadrantA, 'A');
						break;
					case 1:
						processSeqInQuadrant(seq, sLen, quadrantT, 'T');
						break;
					case 2:
						processSeqInQuadrant(seq, sLen, quadrantC, 'C');
						break;
					case 3:
						processSeqInQuadrant(seq, sLen, quadrantG, 'G');
						break;
				}
			}
		}
	}

	// Write result
	if(!writeOutput()){
		return false;
	}
	return true;
}

void KmerChaosGen::processSeqInQuadrant(const string& seq, int sLen, Matrix& quadrant, char base){
	char lBase = tolower(base);
	char uBase = toupper(base);
	for(int kmerStart = 0; kmerStart <= sLen - kmerSize; kmerStart++){
		if(seq[kmerStart] == lBase || seq[kmerStart] == uBase){
			// Kmer's first base is in this quadrant
			// Find where
			int row = 0;
			int col = 0;
			int shift = 1;
			bool invalidBase = false;
			for(int pos = kmerStart + kmerSize - 1; pos > kmerStart; pos--){
				switch(seq[pos]){
					case 'A':
					case 'a':
						break;
					case 'G':
					case 'g':
						col += shift;
						break;
					case 'C':
					case 'c':
						row += shift;
						break;
					case 'T':
					case 't':
						col += shift;
						row += shift;
						break;
					default:
						invalidBase = true;
				}
				shift *= 2;
				if(invalidBase){
					break;
				}
			}
			if(!invalidBase){
				quadrant[row][col] += 1;
			}
		}
	}
}

/*** Finalise results to file
**/
bool KmerChaosGen::writeOutput(){

	if(!outFile.is_open()){
		cerr << "Output file is not open for writing! Can't output results.\n";
		return false;
	}

	for(int row=0; row<quadrantWidth; row++){
		for(int col=0; col<quadrantWidth; col++){
			outFile << quadrantA[row][col] << ",";
		}
		for(int col=0; col<quadrantWidth; col++){
			outFile << quadrantG[col][row] << ",";
		}
		outFile << "\n";
	}
	for(int row=0; row<quadrantWidth; row++){
		for(int col=0; col<quadrantWidth; col++){
			outFile << quadrantC[row][col] << ",";
		}
		for(int col=0; col<quadrantWidth; col++){
			outFile << quadrantT[col][row] << ",";
		}
		outFile << "\n";
	}

	outFile.close();
	ready = false;
	cout << "Done!" << endl;
}

