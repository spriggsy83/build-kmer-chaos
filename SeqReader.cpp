#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "SeqReader.h"
using namespace std;

/* By Andrew Spriggs, CSIRO Ag&Food, 2018 */
/* andrew.spriggs@csiro.au */
/* https://github.com/spriggsy83 */

/*** Construct a SeqReader for (filename) and open file ready for reading. 
**/
SeqReader::SeqReader(const string& aFilename){
	currLen = 0;
	currSeq.clear();
	currID.clear();
	currQual.clear();
	nextID.clear();
	reachedEnd = false;
	
	filename = aFilename;
	bool gzipFile = false;

	if(filename.find("gz", filename.length()-3) != string::npos || 
			filename.find("GZ", filename.length()-3) != string::npos){
		gzipFile = true;
	}
	
	if(gzipFile){
		fileifs.open(filename.c_str(), ios_base::in | ios_base::binary);
	}else{
		fileifs.open(filename.c_str(), ios_base::in );
	}
	if (!fileifs.is_open()){
		cerr << "Unable to open file " << filename << "!\n";
		fileOpen = false;
	}else if (!fileifs.good()){
		cerr << "File " << filename << " is empty!\n";
		fileifs.close();
		fileOpen = false;
	}
	try {
		infile.reset();
		if(gzipFile){
			infile.push(boost::iostreams::gzip_decompressor());
		}
		infile.push(fileifs);
		fileOpen = true;

		char startChar = infile.peek();
		switch (startChar){
			case '@':
				mode = 0; // FASTQ
				break;
			case '>':
				mode = 1; // FASTA
				break;
			default:
				cerr << "File " << filename << " does not appear to be a valid sequence format!\n";
				fileifs.close();
				infile.reset();
				fileOpen = false;
		}
	}
	catch(const boost::iostreams::gzip_error& e) {
		cerr << "Error while reading file " << filename << endl;
		cerr << e.what() << endl;
		fileifs.close();
		infile.reset();
		fileOpen = false;
	}
	return;
}

SeqReader::~SeqReader(){
	currLen = 0;
	if(fileOpen){
		fileifs.close();
		infile.reset();
	}
	fileOpen = false;
}

/*** Fetches the next sequence from file into memory. Returns false if EOF or no sequence read. 
**/
bool SeqReader::nextSeq(){
	
	if(!fileOpen || reachedEnd){
		return false;
	}
	
	currSeq.clear();
	currID.clear();
	currQual.clear();
	currLen = 0;
	
	switch (mode){
		case 0:
			return nextSeqFastq();
		case 1:
		default:
			return nextSeqFasta();
	}
}
	
bool SeqReader::nextSeqFasta(){
	
	long long seqStartMarker;
	bool reachedNext = false;
	bool firstLine = true;
	string lineTmp;
	if(nextID.length() > 0){
		currID = nextID;
		firstLine = false;
	}
	
	try{
		while (infile.good() && !reachedNext){
			getline(infile, lineTmp);
			if(lineTmp.length() > 0){
				if(firstLine){
					if(lineTmp[0] == '>' && lineTmp.length() > 1){
						currID = lineTmp.substr(1);
						firstLine = false;
							// If a space or tab, remove to keep basic ID
						std::size_t spacePos = currID.find_first_of(" \t");
						if(spacePos != std::string::npos){
							currID.replace(spacePos, string::npos, "");
						}
					}else{
						cerr << "File " << filename << " not in valid fasta format!\n";
						return false;
					}
				}else{
					if(lineTmp[0] == '>'){
						nextID = lineTmp.substr(1);
						reachedNext = true;
							// If a space or tab, remove to keep basic ID
						std::size_t spacePos = nextID.find_first_of(" \t");
						if(spacePos != std::string::npos){
							nextID.replace(spacePos, string::npos, "");
						}
					}else{
						currSeq.append(lineTmp);
					}
				}
			}
		}
		if(infile.eof()){
			reachedEnd = true;
		}
		currLen = currSeq.length();
		if(currLen == 0 || currID.length() == 0){
			cerr << "File " << filename << " not in valid fasta format or a sequence was of zero length!\n";
			cerr << currID << "\n" << currSeq << "\n";
			currSeq.clear();
			currID.clear();
			currLen = 0;
			return false;
		}
	}
	catch(const boost::iostreams::gzip_error& e) {
		cerr << "Error while reading .gz file " << filename << endl;
		cerr << e.what() << endl;
		reachedEnd = true;
		return false;
	}
	return true;
}

bool SeqReader::nextSeqFastq(){
	
	long long seqStartMarker;
	bool reachedNext = false;
	string lineTmp;
	int fqLineNum = 0;
	if(nextID.length() > 0){
		currID = nextID;
		fqLineNum++;
	}
	
	try{
		while (infile.good() && !reachedNext){
			getline(infile, lineTmp);
			if(lineTmp.length() > 0){
				if(fqLineNum == 0){
					if(lineTmp[0] == '@' && lineTmp.length() > 1){
						currID = lineTmp.substr(1);
					}else{
						cerr << "File " << filename << " not in valid fastq format!\n";
						cerr << "Invalid line, expecting ID: " << lineTmp << "\n";
						return false;
					}
				}else if(fqLineNum == 1){
					if(isalpha(lineTmp[0])){
						currSeq = lineTmp;
					}else{
						cerr << "File " << filename << " not in valid fastq format!\n";
						cerr << "Invalid line, expecting sequence: " << lineTmp << "\n";
						return false;
					}
				}else if(fqLineNum == 2){
					if(lineTmp[0] != '+'){
						cerr << "File " << filename << " not in valid fastq format!\n";
						cerr << "Invalid line, expecting '+': " << lineTmp << "\n";
						return false;
					}
				}else if(fqLineNum == 3){
					currQual = lineTmp;
				}else if(fqLineNum == 4){
					if(lineTmp[0] == '@' && lineTmp.length() > 1){
						reachedNext = true;
						nextID = lineTmp.substr(1);
					}else{
						cerr << "File " << filename << " not in valid fastq format!\n";
						cerr << "Invalid line, expecting ID: " << lineTmp << "\n";
						return false;
					}
				}
				fqLineNum++;
			}
		}
		if(infile.eof()){
			reachedEnd = true;
		}
		currLen = currSeq.length();
		if(currLen == 0 || currID.length() == 0){
			if(currLen == 0){
				cerr << "File " << filename << " not in valid fastq format; a sequence was of zero length!\n";
			}else{
				cerr << "File " << filename << " not in valid fastq format; missing a sequence ID?!\n";
			}
			currSeq.clear();
			currID.clear();
			currQual.clear();
			currLen = 0;
			return false;
		}
	}
	catch(const boost::iostreams::gzip_error& e) {
		cerr << "Error while reading .gz file " << filename << endl;
		cerr << e.what() << endl;
		reachedEnd = true;
		return false;
	}
	return true;
}

/*** Returns the last sequence string fetched from the file.
**/
string SeqReader::getSeq() const{
	return currSeq;
}

/*** Returns the last sequence ID fetched from the file.
**/
string SeqReader::getSeqID() const{
	return currID;
}

/*** Returns the length of the last sequence fetched from the file.
**/
int SeqReader::getSeqLen() const{
	return currLen;
}

/*** Returns the quality scores of the last sequence fetched from the file (if present).
**/
string SeqReader::getSeqQual() const{
	return currQual;
}

/*** Returns the numeric ID for the file format of the opened file.
*** 0 = FASTQ, 1 = FASTA
**/
int SeqReader::getFileMode() const{
	return mode;
}

/*** Returns the name of the file format of the opened file.
**/
string SeqReader::fileModeString() const{
	string result;
	switch (mode){
		case 0:
			result = "fastq";
			break;
		case 1:
		default:
			result = "fasta";
	}
	return result;
}

/*** Returns full sequence information of the last sequence fetched from the file, as a string.
** Returns sequence ID, sequence and quality scores (if present) in format of file.
**/
string SeqReader::toString() const{
	stringstream result;
	
	switch (mode){
		case 0:
			result << "@" << currID << "\n";
			result << currSeq << "\n";
			result << "+\n" << currQual << "\n";
			break;
		case 1:
		default:
			result << ">" << currID << "\n";
			int printStart = 0;
			int printEnd = 59;
			do{
				if(printEnd >= currLen){
					result << currSeq.substr(printStart) << "\n";
				}else{
					result << currSeq.substr(printStart, printEnd-printStart+1) << "\n";
				}
				printStart += 60;
				printEnd += 60;
			}while(printStart < currLen);
	}
	return result.str();
}


/*** Returns a sub-sequence from the last sequence fetched from the file.
** Start and End are inclusive and count from 0.
**/
string SeqReader::getSubseq(const int ssStart, const int ssEnd) const{
	int start = ssStart;
	int end = ssEnd;
	if(start < 0){
		start = abs(start);
	}
	if(end < 0){
		end = abs(end);
	}
	if(start > end){
		int temp = end;
		end = start;
		start = temp;
	}
	if(start > currLen || start < 1 || end < 1){
		return "";
	}
	return currSeq.substr(start-1, end-start+1);
}

/*** Returns a sub-sequence from the last sequence fetched from the file, with ID, etc., in the format of file.
** Start and End are inclusive and count from 1.
**/
string SeqReader::toStringSubseq(const int ssStart, const int ssEnd) const{
	int start = ssStart;
	int end = ssEnd;
	if(start < 0){
		start = abs(start);
	}
	if(end < 0){
		end = abs(end);
	}
	if(start > end){
		int temp = end;
		end = start;
		start = temp;
	}
	if(start > currLen || start < 1 || end < 1){
		return "";
	}
	
	string newSeq = currSeq.substr(start-1, end-start+1);
	stringstream result;
	
	switch (mode){
		case 0:
			result << "@" << currID << "\n";
			result << newSeq << "\n";
			result << "+\n" << currQual.substr(start-1, end-start+1) << "\n";
			break;
		case 1:
		default:
			result << ">" << currID << ":" << start << "-" << end << "\n";
			int printStart = 0;
			int printEnd = 59;
			do{
				if(printEnd >= newSeq.length()){
					result << newSeq.substr(printStart) << "\n";
				}else{
					result << newSeq.substr(printStart, printEnd-printStart+1) << "\n";
				}
				printStart += 60;
				printEnd += 60;
			}while(printStart < newSeq.length());
	}
	return result.str();
}

/*** Returns the reverse complement of the last sequence string fetched from the file.
**/
string SeqReader::revComp() const{
	string revSeq;
	for(int i=currSeq.length()-1; i>=0; i--){
		switch (currSeq[i]){
			case 'A':
				revSeq.push_back('T');
				break;
			case 'T':
				revSeq.push_back('A');
				break;
			case 'C':
				revSeq.push_back('G');
				break;
			case 'G':
				revSeq.push_back('C');
				break;
			case 'R':
				revSeq.push_back('Y');
				break;
			case 'Y':
				revSeq.push_back('R');
				break;
			case 'M':
				revSeq.push_back('K');
				break;
			case 'K':
				revSeq.push_back('M');
				break;
			case 'S':
				revSeq.push_back('S');
				break;
			case 'W':
				revSeq.push_back('W');
				break;
			case 'B':
				revSeq.push_back('V');
				break;
			case 'D':
				revSeq.push_back('H');
				break;
			case 'H':
				revSeq.push_back('D');
				break;
			case 'V':
				revSeq.push_back('B');
				break;
			case 'a':
				revSeq.push_back('t');
				break;
			case 't':
				revSeq.push_back('a');
				break;
			case 'c':
				revSeq.push_back('g');
				break;
			case 'g':
				revSeq.push_back('c');
				break;
			case 'r':
				revSeq.push_back('y');
				break;
			case 'y':
				revSeq.push_back('r');
				break;
			case 'm':
				revSeq.push_back('k');
				break;
			case 'k':
				revSeq.push_back('m');
				break;
			case 's':
				revSeq.push_back('s');
				break;
			case 'w':
				revSeq.push_back('w');
				break;
			case 'b':
				revSeq.push_back('v');
				break;
			case 'd':
				revSeq.push_back('h');
				break;
			case 'h':
				revSeq.push_back('d');
				break;
			case 'v':
				revSeq.push_back('b');
				break;
			case 'N':
			default:
				revSeq.push_back('N');
		}
	}
	return revSeq;
}

