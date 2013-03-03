#include <map>
#include <set>
#include <vector>
#include <string>
#include <memory>

#include "CountingTrainer.h"
#include "Matrix.h"

using namespace std;

void train_by_counting(HMM& model, vector<string> observations, vector<vector<string>> annotations)
{
    if (model.isFinalized())
        throw invalid_argument("Model must not be finalized!");
    
    Matrix<size_t> A(model.numStates(), model.numStates(), 0);
    vector<size_t> pi(model.numStates(), 0);
    vector<map<string,size_t>> emissions(model.numStates(), map<string,size_t>());
    auto countEmission = [&emissions] (size_t state, string obs) {
        if (emissions[state].count(obs) > 0)
            emissions[state].at(obs)++;
        else
            emissions[state][obs] = 0;
    };
    
    for (int i = 0; i < observations.size(); i++) {
        size_t curpos = 0;
        
        for (int j = 0; j < annotations[i].size()-1; j++) {
            size_t curstate = model.getState(annotations[i][j]);
            
            A(curstate, model.getState(annotations[i][j+1]))++;
            
            string obs = observations[i].substr(curpos, model.stateArity(curstate));
            countEmission(curstate, obs);
            
            curpos += model.stateArity(curstate);
        }
        
        size_t laststate = model.getState(annotations[i][annotations[i].size()-1]);
        string lastobs = observations[i].substr(curpos, model.stateArity(laststate));
        countEmission(laststate, lastobs);
        
        pi[model.getState(annotations[i][0])]++;
    }
    
    auto sumA = [&model, &A] (size_t row) {
        size_t sum = 0;
        for (size_t j = 0; j < model.numStates(); j++)
            sum += A(row, j);
        return sum;
    };
    
    auto sumEmission = [&emissions] (size_t state) {
        size_t sum = 0;
        for (auto emission : emissions[state])
            sum += emission.second;
        return sum;
    };
    
    auto sumPi = [&model, &pi] () {
        size_t sum = 0;
        for (size_t i = 0; i < model.numStates(); i++)
            sum += pi[i];
        return sum;
    };
    
    for (int i = 0; i < model.numStates(); i++) {
        for (int j = 0; j < model.numStates(); j++)
            model.setTransitionProb(i, j, (double) A(i, j) / sumA(i));
        
        for (auto emission : emissions[i])
            model.setEmissionProb(i, emission.first, (double) emission.second / sumEmission(i));
        
        model.setStartProb(i, (double) pi[i] / sumPi());
    }
}
