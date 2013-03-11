#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>

#include "Viterbi.h"
#include "Matrix.h"

using namespace std;

pair<double,vector<size_t>> viterbi(string observation, const HMM& model)
{
    if (!model.isFinalized())
        throw invalid_argument("Model should be finalized!");
    
    // (state, prob)
    Matrix<pair<int,double>> omega(observation.length(), model.numStates(),
                                   make_pair(-1, -numeric_limits<double>::infinity()));

    unordered_map<double, double> logmemory;
    auto ln = [&logmemory] (double arg) {
        if (logmemory.count(arg) > 0)
            return logmemory.at(arg);
        
        double val = log(arg);
        logmemory.insert(make_pair(arg, val));
        return val;
    };
    
    for (size_t i = 0; i < model.numStates(); i++)
        omega(0, i) = make_pair(-1, ln(model.startProb(i)) + ln(model.emissionProb(i, observation.substr(0,1))));
    
    for (size_t l = 1; l < observation.length(); l++) {
        for (size_t i = 0; i < model.numStates(); i++) {
            // Find where we should come from
            pair<int, double> best = make_pair(-1, -numeric_limits<double>::infinity());
            for (auto k : model.incommingStates(i)) {
                if (l < model.stateArity(i))
                    continue;
                    
                double candidate = omega(l - model.stateArity(i) , k).second + ln(model.transitionProb(k, i));
                if (candidate > best.second)
                    best = make_pair(k, candidate);
            }
            
            if (best.first == -1) {
                // State is not possible
                omega(l, i) = make_pair(-1, -numeric_limits<double>::infinity());
            } else {
                // Update current cell with right values
                omega(l, i) = make_pair(best.first,
                                        best.second + ln(model.emissionProb(i, observation.substr(l - model.stateArity(i) + 1, model.stateArity(i)))));
            }
        }
    }
    
    // Final result is now in prev
    pair<int, double> best = make_pair(-1, -numeric_limits<double>::infinity());
    for (int i = 0; i < model.numStates(); i++) {
        double candidate = omega(observation.length()-1, i).second;
        if (candidate > best.second)
            best = make_pair(i, candidate);
    }
    
    if (best.first == -1)
        return make_pair(-numeric_limits<double>::infinity(), vector<size_t>());
    
    // Backtrack
    vector<size_t> stateTrace;
    stateTrace.push_back(best.first);
    size_t pos = observation.length() - 1;
    auto cur = omega(pos, best.first);
    size_t prevState = best.first; // TODO: Could probably be refactored
    while (cur.first != -1) {
        stateTrace.push_back(cur.first);

        pos -= model.stateArity(prevState);
        prevState = cur.first;
        cur = omega(pos, cur.first);
    }
    
    return make_pair(best.second,
                     vector<size_t>(stateTrace.rbegin(), stateTrace.rend()));
}
