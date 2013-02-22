#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <map>

#include "Matrix.h"

using namespace std;

// Natural logarithm function
auto ln = [] (double x) { return log(x); };

class HMM
{
public:
    HMM(vector<string> states, vector<char> symbols)
    : A(Matrix<double>(states.size(), states.size(), 0.)),
    phi(Matrix<double>(states.size(), symbols.size(), 0.)),
    pi(vector<double>(states.size(), 0.)),
    incomming(states.size(), {})
    {
        for (size_t i = 0; i < states.size(); i++)
            stateMap.insert(make_pair(i, states[i]));
        
        for (size_t i = 0; i < symbols.size(); i++)
            symbolMap.insert(make_pair(symbols[i], i));
    }
    
    size_t symbols() const { return symbolMap.size(); }
    size_t states() const { return stateMap.size(); }
    
    /**
     * Apply log to all probabilities and set log-transformed bit to true.
     */
    void logTransform() {
        A.map(ln);
        phi.map(ln);
        for (size_t i = 0; i < pi.size(); i++)
            pi[i] = ln(pi[i]);
        
        // Hack! Wrong place to do this!
        for (size_t i = 0; i < states(); i++) {
            for (size_t j = 0; j < states(); j++) {
                if (transitionProb(i, j) > -numeric_limits<double>::infinity())
                    incomming[j].push_back(i);
            }
        }
        
        logTransformed = true;
    }
    
    /**
     * The probability Pr[symbol | state].
     */
    inline double emissionProb(char symbol, size_t state) const {
        return phi(state, symbolMap.at(symbol));
    }
    
    /**
     * The probability of going from state 'from' to 'to'.
     */
    inline double transitionProb(size_t from, size_t to) const {
        return A(from, to);
    }
    
    vector<size_t> incommingStates(size_t state) const {
        return incomming[state];
    }
    
    /**
     * Returns the name of a given state.
     */
    string stateName(size_t state) const {
        return stateMap.at(state);
    }
    
    bool isLogTransformed() const { return logTransformed; }
    
    Matrix<double> A;
    vector<double> pi;
    Matrix<double> phi;
    
    static unique_ptr<HMM> loadFromStream(ifstream& stream);
    
    void toDot(ofstream& stream);
    
private:
    map<char, size_t> symbolMap;
    map<size_t, string> stateMap;
    
    vector<vector<size_t>> incomming;
    
    bool logTransformed = false;
};
