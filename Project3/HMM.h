#pragma once

#include <vector>
#include <string>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <fstream>

#include "Matrix.h"

using namespace std;

class State
{
public:
    State(string label, size_t d) : label(label), d(d) { }
    
    string getLabel() const { return label; }
    
    size_t emissionArity() const { return d; }
    
    void setEmissionProb(string obs, double prob) {
        if (obs.length() != d)
            throw invalid_argument("Wrong length of observation!");
        
        emissionProbs[obs] = prob;
    }
    
    double getEmissionProb(string obs) const {
        if (obs.length() != d)
            throw invalid_argument("Wrong length of observation!");
        
        if (emissionProbs.count(obs) > 0)
            return emissionProbs.at(obs);
        return 0;
    }
    
    map<string,double> getEmissions() const {
        return emissionProbs;
    }
    
private:
    string label;
    size_t d; // Symbols to emit
    map<string, double> emissionProbs; // Emission probs for a given observation
};

class HMM
{
public:
    HMM(vector<State> states) : states(states),
                                A(Matrix<double>(states.size(), states.size(), 0.)),
                                pi(vector<double>(states.size(), 0.)),
                                finalized(false),
                                incomming(states.size(), {})
    {
        for (int i = 0; i < states.size(); i++) {
            if (stateLabels.count(states[i].getLabel()) > 0)
                throw invalid_argument("Label used twice!");
            stateLabels[states[i].getLabel()] = i;
        }
    }
    
    double transitionProb(size_t from, size_t to) const {
        return A(from, to);
    }
    
    void setTransitionProb(size_t from, size_t to, double prob) {
        if (finalized)
            throw invalid_argument("Model is finalized!");
        
        A(from, to) = prob;
    }
    
    void setTransitionProb(string from, string to, double prob) {
        setTransitionProb(getState(from), getState(to), prob);
    }
    
    double emissionProb(size_t state, string obs) const {
        if (states[state].emissionArity() != obs.length())
            return 0;
        
        return states[state].getEmissionProb(obs);
    }
    
    size_t stateArity(size_t state) const {
        return states[state].emissionArity();
    }
    
    void setEmissionProb(size_t state, string obs, double prob) {
        if (finalized)
            throw invalid_argument("Model is finalized!");
        
        states[state].setEmissionProb(obs, prob);
    }
    
    void setEmissionProb(string state, string obs, double prob) {
        setEmissionProb(getState(state), obs, prob);
    }
    
    void setEmissionProb(string state, vector<string> obs, vector<double> probs) {
        if (obs.size() != probs.size())
            throw invalid_argument("Wrong lengths!");
        
        for (int i = 0; i < obs.size(); i++)
            setEmissionProb(getState(state), obs[i], probs[i]);
    }
    
    double startProb(size_t state) const {
        return pi[state];
    }
    
    void setStartProb(size_t state, double prob) {
        if (finalized)
            throw invalid_argument("Model is finalized!");
        
        pi[state] = prob;
    }
    
    void setStartProb(string state, double prob) {
        setStartProb(getState(state), prob);
    }
    
    vector<size_t> incommingStates(size_t state) const {
        return incomming[state];
    }
    
    size_t numStates() const { return states.size(); }
    
    size_t getState(string label) const {
        if (stateLabels.count(label) > 0)
            return stateLabels.at(label);
        throw invalid_argument("Undefined label!");
    }
    
    string stateLabel(size_t state) const {
        return states[state].getLabel();
    }
    
    void finalize() {
        // TODO: Check validity of all probabilities
        
        for (size_t i = 0; i < numStates(); i++) {
            for (size_t j = 0; j < numStates(); j++) {
                if (transitionProb(i, j) > 0)
                    incomming[j].push_back(i);
            }
        }
        
        finalized = true;
    }
    
    void unlock() {
        finalized = false;
        for (size_t i = 0; i < states.size(); i++)
            incomming[i].clear();
    }
    
    bool isFinalized() const {
        return finalized;
    }
    
    void toDot(ofstream& stream) const {
        if (!finalized)
            throw runtime_error("Model should be finalized!");
        
        stream << "digraph foo {" << endl;
        
        for (int i = 0; i < states.size(); i++) {
            stream << states[i].getLabel() << "[label=\"";
            for (auto emission : states[i].getEmissions()) {
                stream << emission.first << ": " << emission.second << "\\n";
            }
            stream << "\"];" << endl;
            
            for (auto incomming : incommingStates(getState(states[i].getLabel()))) {
                stream << states[incomming].getLabel() << " -> "
                       << states[i].getLabel() << " [label=\"" << transitionProb(incomming, i) << "\"];" << endl;
            }
        }
        
        stream << "}" << endl;
    }
    
private:
    vector<State> states;
    Matrix<double> A;
    vector<double> pi;
    unordered_map<string, size_t> stateLabels;
    
    vector<vector<size_t>> incomming;
    
    bool finalized;
};
