#ifndef SEQREADER_H
#define SEQREADER_H

#include <fstream>
#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

/*** Allows easy parsing of sequence files.
** Currently supports FASTQ and FASTA formats
**/
class SeqReader {
	ifstream fileifs;
	boost::iostreams::filtering_istream infile;
	string filename; //!< Filename of sequence file to be read by this SeqReader
	string currSeq; //!< The last sequence fetched from the file
	string currID; //!< The sequence ID of the last sequence fetched from the file
	string currQual; //!< The quality scores of the last sequence fetched from the file
	int currLen; //!< The sequence length of the last sequence fetched from the file
	bool fileOpen; //!< Is the associated file open and good for reading? true/false
	bool reachedEnd; //!< Has the end of the associated file been reached? true/false
	int mode; //!< The file format of the associated file, as 0 = FASTQ, 1 = FASTA
	string nextID; //!< The sequence ID of the next-to-be-read sequence
	
  public:
	  /*** Construct a SeqReader for (filename) and open file ready for reading. **/
	SeqReader(const string&);
	~SeqReader();
		/*** Fetches the next sequence from file into memory. Returns false if EOF or no sequence read. **/
	bool nextSeq();
		/*** Returns the last sequence string fetched from the file.
		**/
	string getSeq() const;
		/*** Returns the last sequence ID fetched from the file.
		**/
	string getSeqID() const;
		/*** Returns the length of the last sequence fetched from the file. **/
	int getSeqLen() const;
		/*** Returns the quality scores of the last sequence fetched from the file (if present). **/
	string getSeqQual() const;
		/*** Returns the numeric ID for the file format of the opened file.
		*** 0 = FASTQ, 1 = FASTA **/
	int getFileMode() const;
		/*** Returns the name of the file format of the opened file. **/
	string fileModeString() const;
		/*** Returns full sequence information of the last sequence fetched from the file, as a string.
		** Returns sequence ID, sequence and quality scores (if present) in format of file. **/
	string toString() const;
		/*** Returns a sub-sequence from the last sequence fetched from the file.
		** Start and End are inclusive and count from 1. **/
	string getSubseq(const int start, const int end) const;
		/*** Returns a sub-sequence from the last sequence fetched from the file, with ID, etc., in the format of file.
		** Start and End are inclusive and count from 1. **/
	string toStringSubseq(const int start, const int end) const;
		/*** Returns the reverse complement of the last sequence string fetched from the file.
		**/
	string revComp() const;
  private:
  	bool nextSeqFastq();
  	bool nextSeqFasta();
};

#endif
