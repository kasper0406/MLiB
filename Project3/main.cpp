#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <random>
#include <functional>
#include <memory>

#include "HMM.h"
#include "Fasta.h"
#include "Viterbi.h"
#include "CountingTrainer.h"
#include "ViterbiTrainer.h"
#include "EMTrainer.h"
#include "SimpleParser.h"
#include "ForwardBackward.h"

using namespace std;

template<class Parser>
vector<string> parse_observation(string observation, string annotation)
{
    Parser parser;
    vector<string> states;
    
    auto add = [&states] (vector<string> visitedStates) {
        for (auto state : visitedStates)
            states.push_back(state);
    };
    auto isReverse = [&annotation] (int i) { return toupper(annotation[i]) == 'R'; };
    char startDirection;
    
    for (int i = 0; i < observation.length();) {
        if (annotation[i] == 'N') {
            states.push_back( parser.handleNoncoding(observation[i]) );
            startDirection = 'N';
            i++;
            continue;
        }
        
        startDirection = toupper(annotation[i]);
        
        // Read start codon
        string startCodon = observation.substr(i, 3);
        if (isReverse(i))
            add( parser.handleReverseStart(startCodon) );
        else
            add( parser.handleStart(startCodon) );
        
        for (i += 3; i+6 < observation.length() && toupper(annotation[i+3]) == startDirection; i += 3) {
            // Read intermediate codons
            string codon = observation.substr(i, 3);
            if (isReverse(i))
                add( parser.handleReverseCoding(codon) );
            else
                add( parser.handleCoding(codon) );
        }
        
        // Read end codon
        string endCodon = observation.substr(i, 3);
        if (!endCodon.empty()) {
            if (isReverse(i))
                add( parser.handleReverseEnd(endCodon) );
            else
                add( parser.handleEnd(endCodon) );
        }
        
        i += 3;
    }
    
    return states;
}

vector<double> random_distr(unsigned int size, unsigned int seed)
{
    mt19937 generator(seed);
    uniform_int_distribution<int> distribution(1, 100);
    auto random = bind(distribution, generator);
    
    vector<unsigned int> numbers(size, 0);
    unsigned int sum = 0;
    for (int i = 0; i < size; i++) {
        numbers[i] = random();
        sum += numbers[i];
    }
    
    vector<double> distr(size, 0.);
    for (int i = 0; i < size; i++)
        distr[i] = (double) numbers[i] / sum;
    return distr;
}

HMM build_model() {
    vector<State> states;
    states.push_back(State("N", 1));
    states.push_back(State("S", 3));
    states.push_back(State("E", 3));
    states.push_back(State("C", 3));
    states.push_back(State("RS", 3));
    states.push_back(State("RE", 3));
    states.push_back(State("RC", 3));
    
    return HMM(states);
}

HMM build_model_with_transitions() {
    HMM model = build_model();
    
    model.setEmissionProb("N", {"A","C","G","T"}, {.25,.25,.25,.25});
    for (auto state : {"S","E","C","RS","RE","RC"}) {
        char symbols[] = { 'A', 'C', 'G', 'T' };
        auto distr = random_distr(4 * 4 * 4, 0xDEADBEEF + (unsigned int)hash<string>()(string(state)));
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    stringstream obs;
                    obs << symbols[i] << symbols[j] << symbols[k];
                    model.setEmissionProb(state, obs.str(), distr[4*4*i+4*j+k]);
                }
            }
        }
    }
    
    model.setTransitionProb("N", "N", .9);
    
    model.setTransitionProb("N", "S", .05);
    model.setTransitionProb("S", "C", 1);
    model.setTransitionProb("C", "C", .95);
    model.setTransitionProb("C", "E", .05);
    model.setTransitionProb("E", "N", 1);
    
    model.setTransitionProb("N", "RS", .05);
    model.setTransitionProb("RS", "RC", 1);
    model.setTransitionProb("RC", "RC", .95);
    model.setTransitionProb("RC", "RE", .05);
    model.setTransitionProb("RE", "N", 1);
    
    model.setStartProb("N", 1);
    
    return model;
}

HMM test_model() {
    /*
    vector<State> states;
    for (int i = 1; i <= 7; i++) {
        stringstream ss;
        ss << i;
        states.push_back(State(ss.str(), 1));
    }
    
    HMM hmm(states);
    hmm.setTransitionProb("4", "4", 0.9);
    hmm.setTransitionProb("4", "3", 0.05);
    hmm.setTransitionProb("4", "5", 0.05);
    hmm.setTransitionProb("3", "2", 1);
    hmm.setTransitionProb("2", "1", 1);
    hmm.setTransitionProb("1", "3", 0.9);
    hmm.setTransitionProb("1", "4", 0.1);
    hmm.setTransitionProb("5", "6", 1);
    hmm.setTransitionProb("6", "7", 1);
    hmm.setTransitionProb("7", "5", 0.9);
    hmm.setTransitionProb("7", "4", 0.1);
    
    hmm.setEmissionProb("1", {"A","C","G","T"}, {0.3,0.25,0.25,0.20});
    hmm.setEmissionProb("2", {"A","C","G","T"}, {0.2,0.35,0.15,0.30});
    hmm.setEmissionProb("3", {"A","C","G","T"}, {0.4,0.15,0.2,0.25});
    hmm.setEmissionProb("4", {"A","C","G","T"}, {0.25,0.25,0.25,0.25});
    hmm.setEmissionProb("5", {"A","C","G","T"}, {0.2,0.4,0.3,0.1});
    hmm.setEmissionProb("6", {"A","C","G","T"}, {0.3,0.2,0.3,0.2});
    hmm.setEmissionProb("7", {"A","C","G","T"}, {0.15,0.30,0.20,0.35});
    
    hmm.setStartProb("4", 1);
     */
    
    vector<State> states;
    states.push_back(State("NC", 1));
    states.push_back(State("C", 3));
    states.push_back(State("R", 3));
    HMM hmm(states);
    
    char symbols[] = { 'A', 'C', 'G', 'T' };
    double prob1[] = {0.3,0.25,0.25,0.20};
    double prob2[] = {0.2,0.35,0.15,0.30};
    double prob3[] = {0.4,0.15,0.2,0.25};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                stringstream obs;
                obs << symbols[i] << symbols[j] << symbols[k];
                hmm.setEmissionProb("C", obs.str(), prob3[i] * prob2[j] * prob1[k]);
            }
        }
    }
    
    hmm.setEmissionProb("NC", {"A","C","G","T"}, {0.25,0.25,0.25,0.25});
    
    double prob5[] = {0.2,0.4,0.3,0.1};
    double prob6[] = {0.3,0.2,0.3,0.2};
    double prob7[] = {0.15,0.30,0.20,0.35};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                stringstream obs;
                obs << symbols[i] << symbols[j] << symbols[k];
                hmm.setEmissionProb("R", obs.str(), prob5[i] * prob6[j] * prob7[k]);
            }
        }
    }
    
    hmm.setStartProb("NC", 1);
    
    hmm.setTransitionProb("NC", "NC", 0.9);
    hmm.setTransitionProb("NC", "C", 0.05);
    hmm.setTransitionProb("NC", "R", 0.05);
    hmm.setTransitionProb("C", "C", 0.9);
    hmm.setTransitionProb("C", "NC", 0.1);
    hmm.setTransitionProb("R", "R", 0.9);
    hmm.setTransitionProb("R", "NC", 0.1);
    
    return hmm;
}

int main(int argc, const char * argv[])
{
    cout << "Loading files..." << endl;
    auto observations = read_seqs_from_files({"genome1.fa","genome2.fa","genome3.fa","genome4.fa","genome5.fa"});
    
    auto annotations = read_seqs_from_files({"annotation1.fa","annotation2.fa","annotation3.fa","annotation4.fa","annotation5.fa"});
    cout << "Parsing observations..." << endl;
    vector<vector<string>> parsed;
    for (int i = 0; i < observations.size(); i++)
        parsed.push_back( parse_observation<SimpleParser>(observations[i], annotations[i]) );
    
    cout << "Building and traning model..." << endl;
    
    HMM model = build_model();
    train_by_counting(model, observations, parsed);
    
    cout << "Stop!" << endl;
    
    // HMM model = build_model_with_transitions();
    //train_by_viterbi(model, observations, 1);

    /*
    ifstream input("predictions/model_bw_40.dot", ifstream::in);
    HMM model = HMM::loadFromDot(input);
    model.setStartProb(0, 1);
    input.close();
     */
    
    model.finalize();
    
    auto toBePredicted = read_seqs_from_files({"genome6.fa","genome7.fa","genome8.fa","genome9.fa","genome10.fa","genome11.fa"});
    
    const int iterations = 20;
    
    /*
    for (int i = 1; i <= iterations; i++) {
        cout << "Viterbi " << i << endl;
        train_by_viterbi(model, observations, 1);
    }
     */
    
    for (int i = 1; i <= iterations; i++) {
        model.unlock();
        
        cout << "Running iteration " << i << " of viterbi training." << endl;
        train_by_baumwelch(model, observations);
        // train_by_viterbi(model, observations, 1);
    
        model.finalize();
    
        cout << "Writing model to dot file..." << endl;
        stringstream modelname;
        modelname << "predictions/model_bwvit_" << i << ".dot";
        ofstream out(modelname.str(), ofstream::out);
        model.toDot(out);
        out.close();
    
        cout << "Running Viterbi..." << endl;
        for (int j = 0; j < toBePredicted.size(); j++) {
            double probability;
            vector<size_t> trace;
            
            tie(probability, trace) = viterbi(toBePredicted[j], model);
            
            cout << "Writing trace..." << endl;
            
            stringstream file;
            file << "predictions/bwvit" << i << "_" << (j+6) << ".fa";
            
            ofstream outpred(file.str(), ofstream::out);
            for (auto s : trace) {
                switch (model.stateLabel(s)[0]) {
                    case 'N': outpred << 'N'; break;
                    case 'R': outpred << string(3, 'R'); break;
                    default: outpred << string(3, 'C'); break;
                }
            }
            outpred.close();
        }
    }
    
    
    /*
    HMM model = test_model();
    model.finalize();
    
    cout << "Writing model to dot file..." << endl;
    ofstream out("testmodel.dot", ofstream::out);
    model.toDot(out);
    out.close();
    
    double probability;
    vector<size_t> trace;
    tie(probability, trace) = viterbi("GTTTCCCAGTGTATATCGAGGGATACTACGTGCATAGTAACATCGGCCAA", model);
    
    cout << "Prob: " << probability << endl;
    for (auto s : trace) {
        if (s == 0) cout << "4";
        else if (s == 1) cout << "321";
        else if (s == 2) cout << "567";
    }
    cout << endl;
     */
    
    /*
    string observation = "GTTTCCCAGTGTATATCGAGGGATACTACGTGCATAGTAACATCGGCCAA";
    auto res = forward_backward(observation, model);
    double C = 1;
    for (size_t i = 0; i < observation.length(); i++) {
        C *= get<0>(res)[i];
        for (size_t state = 0; state < model.numStates(); state++) {
            cout << get<1>(res)(i, state) * C << "\t";
        }
        cout << endl;
    }
    
    cout << endl << endl;
    
    C = 1;
    for (long i = observation.length() - 1; i >= 0; i--) {
        if (i != observation.length() - 1)
            C *= get<0>(res)[i+1];
        for (size_t state = 0; state < model.numStates(); state++) {
            cout << get<2>(res)(i, state) * C << "\t";
        }
        cout << endl;
    }
     */
    
    /*
    cout << "Running Baum-Welch training." << endl;
    train_by_baumwelch(model, {observation});
    
    model.finalize();
    
    cout << "Writing model to dot file..." << endl;
    ofstream outbw("testmodel_bw.dot", ofstream::out);
    model.toDot(outbw);
    outbw.close();
     */
    
    return 0;
}
