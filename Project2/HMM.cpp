#include <fstream>
#include <vector>
#include <memory>

#include "HMM.h"

using namespace std;

unique_ptr<HMM> HMM::loadFromStream(ifstream &file)
{
    // Read states
    string tmp;
    file >> tmp;
    
    size_t stateslen;
    file >> stateslen;
    
    vector<string> states;
    for (size_t i = 0; i < stateslen; i++) {
        string state;
        file >> state;
        states.push_back(state);
    }
    
    // Read observables
    file >> tmp;
    
    size_t observableslen;
    file >> observableslen;
    
    vector<char> observables;
    for (size_t i = 0; i < observableslen; i++) {
        char observation;
        file >> observation;
        observables.push_back(observation);
    }
    
    unique_ptr<HMM> model(new HMM(states, observables));
    
    // Initial probabilities
    file >> tmp;
    for (int i = 0; i < stateslen; i++)
        file >> model->pi[i];
    
    // Transition probabilities
    file >> tmp;
    for (int i = 0; i < stateslen; i++) {
        for (int j = 0; j < stateslen; j++)
            file >> model->A(i, j);
    }
    
    // Emission probabilities
    file >> tmp;
    for (int i = 0; i < stateslen; i++) {
        for (int j = 0; j < observableslen; j++)
            file >> model->phi(i, j);
    }
    
    model->logTransform();
    
    return model;
}

void HMM::toDot(ofstream &stream)
{
    stream << "digraph foo {" << endl;
    
    for (auto state : stateMap) {
        stream << state.second << "[label=\"";
        for (auto symbol : symbolMap) {
            stream << symbol.first << ": " << emissionProb(symbol.first, state.first) << "\\n";
        }
        stream << "\"];" << endl;
        
        for (auto incomming : incommingStates(state.first)) {
            stream << stateMap.at(incomming) << " -> " << state.second << " [label=\"" << transitionProb(incomming, state.first) << "\"];" << endl;
        }
    }
    
    stream << "}" << endl;
}
