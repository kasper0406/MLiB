#ifndef __Project2__fasta__
#define __Project2__fasta__

#include <vector>
#include <fstream>
#include <string>

using namespace std;

vector<pair<string,string>> read_fasta_from_stream(ifstream& stream);
vector<string> read_seqs_from_files(vector<string> files);

#endif /* defined(__Project2__fasta__) */
