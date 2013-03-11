#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include "ForwardBackward.h"
#include "Matrix.h"
#include "HMM.h"

using namespace std;

tuple<vector<double>, Matrix<double>, Matrix<double>> forward_backward(string obs, const HMM& model)
{
    if (!model.isFinalized())
        throw runtime_error("Model should be finalized!");
    
    // Forward algorithm
    Matrix<double> forward(obs.length(), model.numStates(), 0);
    vector<double> cs(obs.length(), 0);
    
    // Calculate c1
    for (size_t state = 0; state < model.numStates(); state++)
        cs[0] += model.startProb(state) * model.emissionProb(state, obs.substr(0,1));
    // Base case
    for (size_t state = 0; state < model.numStates(); state++)
        forward(0, state) = model.startProb(state) * model.emissionProb(state, obs.substr(0,1)) / cs[0];

    // Recursion
    for (size_t i = 1; i < obs.length(); i++) {
        vector<double> delta(model.numStates(), 0);
        for (size_t state = 0; state < model.numStates(); state++) {
            if (i < model.stateArity(state))
                continue;
            
            for (auto prevState : model.incommingStates(state)) {
                double val = forward(i - model.stateArity(state), prevState) * model.transitionProb(prevState, state);
                for (size_t k = 1; k < model.stateArity(state); k++)
                    val /= cs[i - k];
                
                delta[state] += val;
            }
            delta[state] *= model.emissionProb(state, obs.substr(i - model.stateArity(state) + 1, model.stateArity(state)));
            
            cs[i] += delta[state];
        }
        
        for (size_t state = 0; state < model.numStates(); state++) {
            forward(i, state) = delta[state] / cs[i];
        }
    }
    
    // Backward algorithm
    Matrix<double> backward(obs.length(), model.numStates(), 0);
    const size_t N = obs.length() - 1;
    for (size_t state = 0; state < model.numStates(); state++)
        backward(N, state) = 1;
    
    for (long i = N - 1; i >= 0; i--) {
        for (size_t state = 0; state < model.numStates(); state++) {
            double prob = 0;
            
            for (auto nextState : model.outgoingStates(state)) {
                if (i + model.stateArity(nextState) > N)
                    continue;
                
                double val = backward(i + model.stateArity(nextState), nextState) * model.transitionProb(state, nextState)
                               * model.emissionProb(nextState, obs.substr(i + 1, model.stateArity(nextState)));
                
                for (size_t k = 0; k < model.stateArity(nextState); k++)
                    val /= cs[i + 1 + k];
                
                prob += val;
            }
            backward(i, state) = prob;
        }
    }
    
    return make_tuple(cs, forward, backward);
}
