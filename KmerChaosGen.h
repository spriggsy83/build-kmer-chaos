#ifndef KMERCHAOS_H
#define KMERCHAOS_H

#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctype.h>
#include <math.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2019 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

typedef vector< unsigned int > MatrixRow;
typedef vector< MatrixRow > Matrix; //!< to hold counts per kmer in a matrix

class KmerChaosGen {
  private:
  	int kmerSize; //!< Kmer size to process from sequence
  	string inFileName;
  	string outFileName;
  	ofstream outFile;
  	bool ready;  //!< Indicates that class has been initialised and processing is ready to run
	
	Matrix quadrantA; //!< counts per kmer in a matrix, A quadrant
	Matrix quadrantT; //!< counts per kmer in a matrix, T quadrant
	Matrix quadrantC; //!< counts per kmer in a matrix, C quadrant
	Matrix quadrantG; //!< counts per kmer in a matrix, G quadrant
	unsigned int quadrantWidth; // The width and height of the result matrix
	
		/*** Actual constructor code, to be called by any constructor forms **/
	void prepareKmerChaosGen(const string& aInFileName, const string& aOutFileName, const int akmerSize);
		/*** Fill in one quadrants worth of kmers **/
	void processSeqInQuadrant(const string& seq, int sLen, Matrix& quadrant, char base);
		/*** Finalise results to file **/
	bool writeOutput();

  public:
		/** Constructor **/
	KmerChaosGen(const string& aInFileName, const string& aOutFileName, const int akmerSize );
	~KmerChaosGen();
	
		/** Launch the full kmer analysis process **/
	bool GenKmerChaos();
};

#endif
