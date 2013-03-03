#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <memory>
#include <fstream>
#include <cassert>

#include "TreeParser.h"
#include "fasta.h"
#include "CountingTrainer.h"
#include "Viterbi.h"

using namespace std;

template<class Parser>
vector<string> parse_observation(string observation, string annotation)
{
    vector<string> states;
    
    auto add = [&states] (tuple<string,string,string> s) {
        states.push_back(get<0>(s));
        states.push_back(get<1>(s));
        states.push_back(get<2>(s));
    };
    auto isReverse = [&annotation] (int i) { return toupper(annotation[i]) == 'R'; };
    char startDirection;

    for (int i = 0; i < observation.length();) {
        if (annotation[i] == 'N') {
            states.push_back( Parser::handleNoncoding(observation[i]) );
            startDirection = 'N';
            i++;
            continue;
        }
        
        startDirection = toupper(annotation[i]);
        
        // Read start codon
        string startCodon = observation.substr(i, 3);
        if (isReverse(i))
            add( Parser::handleReverseStart(startCodon) );
        else
            add( Parser::handleStart(startCodon) );
        
        for (i += 3; i+6 < observation.length() && toupper(annotation[i+3]) == startDirection; i += 3) {
            // Read intermediate codons
            string codon = observation.substr(i, 3);
            if (isReverse(i))
                add( Parser::handleReverseCoding(codon) );
            else
                add( Parser::handleCoding(codon) );
        }
        
        // Read end codon
        string endCodon = observation.substr(i, 3);
        if (isReverse(i))
            add( Parser::handleReverseEnd(endCodon) );
        else
            add( Parser::handleEnd(endCodon) );
        
        i += 3;
    }

    return states;
}

int main(int argc, const char * argv[])
{
    cout << "Loading files..." << endl;
    
    auto observations = read_seqs_from_files({"genome1.fa","genome2.fa","genome3.fa","genome4.fa","genome5.fa"});
    auto annotations = read_seqs_from_files({"annotation1.fa","annotation2.fa","annotation3.fa","annotation4.fa","annotation5.fa"});
    
    cout << "Parsing observations..." << endl;
    
    vector<vector<string>> parsed;
    for (int i = 0; i < observations.size(); i++)
        parsed.push_back( parse_observation<TreeParser>(observations[i], annotations[i]) );
    
    cout << "Building and traning model..." << endl;
    
    unique_ptr<HMM> model = train_by_counting(observations, parsed);
    
    cout << "Writing model to dot file..." << endl;
    
    ofstream out("model.dot", ofstream::out);
    model->toDot(out);
    out.close();
    
    cout << "Running Viterbi..." << endl;
    
    double probability;
    vector<size_t> trace;
    tie(probability, trace) = viterbi(observations[0], *model);
    
    cout << "Writing trace..." << endl;
    ofstream outpred("prediction_foobar_1.fa", ofstream::out);
    for (auto s : trace) {
        switch (model->stateName(s)[0]) {
            case 'N': outpred << 'N'; break;
            case 'R': outpred << 'R'; break;
            default: outpred << 'C'; break;
        }
    }
    outpred.close();
    
    return 0;
}

