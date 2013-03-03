#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>

#include "Fasta.h"

using namespace std;

vector<pair<string,string>> read_fasta_from_stream(ifstream& stream)
{
    vector<pair<string,string>> seqs;
    string seq, line, name;
    while (getline(stream, line)) {
        if (!line.empty() && line[0] == '>') {
            if (!seq.empty()) {
                seqs.push_back(make_pair(name, seq));
                seq = "";
            }
            name = line.substr(1);
        } else if (!line.empty() && line[0] == ';') {
            // Ignore the line...
        } else {
            for (unsigned int i = 0; i < line.length(); i++) {
                if (line[i] != ' ')
                    seq += line[i];
            }
        }
    }
    seqs.push_back(make_pair(name, seq));
    return seqs;
}

vector<string> read_seqs_from_files(vector<string> files)
{
    vector<string> seqs;
    for (string file : files) {
        ifstream stream(file);
        if (!stream.is_open())
            throw runtime_error("Could not find " + file);
        
        auto seqsInFile = read_fasta_from_stream(stream);
        for (auto seqdesc : seqsInFile)
            seqs.push_back(seqdesc.second);
        
        stream.close();
    }
    return seqs;
}
