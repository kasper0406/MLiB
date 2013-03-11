#pragma once

#include <vector>
#include <string>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <fstream>
#include <cmath>
#include <sstream>

#include "Matrix.h"

using namespace std;

class State
{
public:
    State(string label, size_t d) : label(label), d(d) {
        if (d < D)
            emissionProbsVec = vector<double>(pow(4, d), 0);
    }
    
    string getLabel() const { return label; }
    
    size_t emissionArity() const { return d; }
    
    void setEmissionProb(string obs, double prob) {
        if (obs.length() != d)
            throw invalid_argument("Wrong length of observation!");
        
        if (d < D)
            emissionProbsVec[getIndex(obs)] = prob;
        
        if (prob == 0) {
            if (emissionProbs.count(obs) > 0)
                emissionProbs.erase(obs);
        } else
            emissionProbs[obs] = prob;
    }
    
    double getEmissionProb(string obs) const {
        if (obs.length() != d)
            throw invalid_argument("Wrong length of observation!");
        
        if (d < D)
            return emissionProbsVec[getIndex(obs)];
        
        if (emissionProbs.count(obs) > 0)
            return emissionProbs.at(obs);
        return 0;
    }
    
    void resetEmissions() {
        emissionProbs.clear();
        if (d < D) {
            for (int i = 0; i < pow(4, d); i++)
                emissionProbsVec[i] = 0;
        }
    }
    
    unordered_map<string,double> getEmissions() const {
        return emissionProbs;
    }
    
private:
    string label;
    size_t d; // Symbols to emit
    static const size_t D = 5;
    
    unordered_map<string, double> emissionProbs; // Emission probs for a given observation
    vector<double> emissionProbsVec; // Emission probs stored in a vector
    
    inline unsigned int getIndex(char c) const {
        switch (c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default:
                throw runtime_error("Invalid symbol!");
        }
    }
    
    inline size_t getIndex(string obs) const {
        size_t index = 0;
        for (int i = 0; i < obs.length(); i++) {
            index += getIndex(obs[i]) * (1 << (2*i));
        }
        return index;
    }
};

class HMM
{
public:
    HMM(vector<State> states) : states(states),
                                A(Matrix<double>(states.size(), states.size(), 0.)),
                                pi(vector<double>(states.size(), 0.)),
                                finalized(false),
                                incomming(states.size(), {}),
                                outgoing(states.size(), {})
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
    
    vector<size_t> outgoingStates(size_t state) const {
        return outgoing[state];
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
        if (finalized)
            throw runtime_error("Model already finalized!");
        
        for (size_t i = 0; i < numStates(); i++) {
            double emissionProb = 0;
            for (auto emission : states[i].getEmissions())
                emissionProb += emission.second;
            if (abs(emissionProb - 1) > EPSILON) {
                stringstream ss;
                ss << "Emission probs does not sum to 1 in state '" << stateLabel(i) << "'!";
                throw runtime_error(ss.str());
            }
            
            double transProb = 0;
            for (size_t j = 0; j < numStates(); j++)
                transProb += transitionProb(i, j);
            if (abs(transProb - 1) > EPSILON) {
                stringstream ss;
                ss << "Transition probs does not sum to 1 in state '" << stateLabel(i) << "'!";
                throw runtime_error(ss.str());
            }
        }
        
        for (size_t i = 0; i < numStates(); i++) {
            for (size_t j = 0; j < numStates(); j++) {
                if (transitionProb(i, j) > 0) {
                    incomming[j].push_back(i);
                    outgoing[i].push_back(j);
                }
            }
        }
        
        finalized = true;
    }
    
    void unlock() {
        finalized = false;
        for (size_t i = 0; i < states.size(); i++) {
            incomming[i].clear();
            outgoing[i].clear();
        }
    }
    
    void reset() {
        if (finalized)
            throw runtime_error("Model must not be finalized!");
        
        A.map([] (double x) { return 0; });
        for (size_t i = 0; i < numStates(); i++) {
            states[i].resetEmissions();
            pi[i] = 0;
        }
    }
    
    bool isFinalized() const {
        return finalized;
    }
    
    unordered_map<string, double> getEmissions(size_t state) const {
        return states[state].getEmissions();
    }
    
    void toDot(ofstream& stream) const {
        if (!finalized)
            throw runtime_error("Model should be finalized!");
        
        stream << "digraph foo {" << endl;
        
        for (int i = 0; i < states.size(); i++) {
            stream << states[i].getLabel() << "[label=\"";
            for (auto emission : getEmissions(i)) {
                stream << emission.first << ": " << emission.second << "\\n";
            }
            stream << "\",shape=box];" << endl;
            
            for (auto incomming : incommingStates(i)) {
                stream << states[incomming].getLabel() << " -> "
                       << states[i].getLabel() << " [label=\"" << transitionProb(incomming, i) << "\"];" << endl;
            }
        }
        
        stream << "}" << endl;
    }
    
    static HMM loadFromDot(ifstream& is) {
        map<size_t, State> states;
        map<pair<size_t, size_t>, double> transitions;
        vector<map<string, double>> emissions;
        
        map<string, size_t> labelToIndex;
        auto stateIndex = [&labelToIndex, &emissions] (string label) {
            if (labelToIndex.count(label) > 0)
                return labelToIndex.at(label);
            else {
                labelToIndex[label] = labelToIndex.size();
                emissions.push_back(map<string,double>());
                return labelToIndex[label];
            }
        };
        
        string line;
        while (getline(is, line)) {
            if (line.find(" -> ") != string::npos) {
                // Transition line
                size_t arrowStart = line.find(" -> ");
                string from = line.substr(0, arrowStart);
                size_t stop = line.find(" [label=\"");
                string to = line.substr(arrowStart + 4, stop - (arrowStart + 4));
                size_t probstop = line.find("\"];");
                istringstream iss(line.substr(stop + 9, probstop - (stop+9)));
                double prob;
                iss >> prob;
                
                transitions.insert(make_pair(make_pair(stateIndex(from), stateIndex(to)), prob));
            } else if (line.find("[label=\"") != string::npos) {
                // State line
                size_t probsStart = line.find("[label=\"");
                string state = line.substr(0, probsStart);
                size_t labelStop = line.find("\\n\",shape=box];");
                string label = line.substr(probsStart + 8, labelStop - (probsStart + 8));
                size_t d;
                
                size_t pos = 0;
                do {
                    size_t newpos = label.find("\\n", pos);
                    string emission;
                    if (newpos == string::npos)
                        emission = label.substr(pos, newpos);
                    else
                        emission = label.substr(pos, newpos - pos);
                    
                    size_t sep = emission.find(": ");
                    string obs = emission.substr(0, sep);
                    istringstream probstream(emission.substr(sep+2));
                    double prob = 0;
                    probstream >> prob;
                    
                    d = obs.length();
                    emissions.at(stateIndex(state)).insert(make_pair(obs, prob));
                    
                    if (newpos != string::npos)
                        pos = newpos + 2;
                    else
                        pos = newpos;
                } while (pos != string::npos);
                
                states.insert(make_pair(stateIndex(state), State(state, d)));
            } else {
                if (line == "digraph foo {" || line == "}") continue;
                
                throw runtime_error("Illegal line!");
            }
        }
        
        vector<State> orderedStates;
        for (auto s : states)
            orderedStates.push_back(s.second);
        
        HMM model(orderedStates);
        for (auto transition : transitions) {   
            model.setTransitionProb(transition.first.first, transition.first.second, transition.second);
        }
        
        for (int i = 0; i < model.numStates(); i++) {
            for (auto emission : emissions[i])
                model.setEmissionProb(i, emission.first, emission.second);
        }
        
        return model;
    }
    
private:
    vector<State> states;
    Matrix<double> A;
    vector<double> pi;
    unordered_map<string, size_t> stateLabels;
    
    vector<vector<size_t>> incomming, outgoing;
    
    bool finalized;
    
    constexpr static const double EPSILON = 0.000001;
};
