#include <vector>
#include <string>
#include <iostream>

#include "EMTrainer.h"
#include "HMM.h"
#include "ForwardBackward.h"

using namespace std;

void train_by_baumwelch(HMM& model, vector<string> observations)
{
    Matrix<double> A(model.numStates(), model.numStates(), 0);
    vector<double> pi(model.numStates(), 0);

    vector<map<string,double>> emissions(model.numStates(), map<string,double>());
    vector<double> throughStateProbs(model.numStates(), 0);
    
    model.finalize();
    
    for (auto observation : observations) {
        auto FBtable = forward_backward(observation, model);
        
        auto gamma = [&model, &FBtable] (size_t n, size_t state) {
            size_t pos = n + model.stateArity(state) - 1;
            return get<1>(FBtable)(pos, state) * get<2>(FBtable)(pos, state);
        };
        
        for (size_t k = 0; k < model.numStates(); k++) {
            for (size_t n = 0; n < observation.length(); n++) {
                if (n + model.stateArity(k) >= observation.length())
                    continue;
                
                string obs = observation.substr(n, model.stateArity(k));
                double emissionProb = model.emissionProb(k, obs);
                
                if (n > 0) {
                    // Transition probabilities
                    double C = 1;
                    for (int i = 0; i < model.stateArity(k); i++)
                        C *= get<0>(FBtable)[n+i];
                    
                    for (size_t j = 0; j < model.numStates(); j++) {
                        A(j, k) += get<1>(FBtable)(n-1, j) * get<2>(FBtable)(n + model.stateArity(k) - 1, k)
                        * (emissionProb * model.transitionProb(j, k)) / C;
                    }
                }
                
                // Emission probabilities
                double gamma_nk = gamma(n, k);                
                if (emissions[k].count(obs) > 0)
                    emissions[k].at(obs) += gamma_nk;
                else
                    emissions[k][obs] = gamma_nk;
                throughStateProbs[k] += gamma_nk;
            }
            
            pi[k] += gamma(0, k);
        }
    }
    
    model.unlock();
    model.reset();
    
    // Update the model
    auto sumA = [&model, &A] (size_t row) {
        double sum = 0;
        for (size_t j = 0; j < model.numStates(); j++)
            sum += A(row, j);
        return sum;
    };
    
    auto sumPi = [&model, &pi] () {
        double sum = 0;
        for (size_t i = 0; i < model.numStates(); i++)
            sum += pi[i];
        return sum;
    };
    
    for (int i = 0; i < model.numStates(); i++) {
        for (int j = 0; j < model.numStates(); j++)
            model.setTransitionProb(i, j, A(i, j) / sumA(i));
        
        for (auto emission : emissions[i]) {
            model.setEmissionProb(i, emission.first, emission.second / throughStateProbs[i]);
        }
        
        model.setStartProb(i, pi[i] / sumPi());
    }
}
