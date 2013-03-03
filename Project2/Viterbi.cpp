#include <vector>
#include <memory>
#include <string>
#include <stdexcept>

#include "Viterbi.h"
#include "HMM.h"

using namespace std;

class List
{
public:
    List(size_t state, shared_ptr<List> tail) : state(state), prev(tail)
    { }
    
    //private:
    size_t state;
    shared_ptr<List> prev = shared_ptr<List>(nullptr);
};

class Cell
{
public:
    double prob = -numeric_limits<double>::infinity(); // Sort of abuse of probabilites!
    shared_ptr<List> list = shared_ptr<List>(nullptr);
};
    
pair<double,vector<size_t>> viterbi(string observation, const HMM& model)
{
    if (!model.isLogTransformed())
        throw runtime_error("Model should be transformed!");
    
    vector<Cell> prev(vector<Cell>(model.states(), Cell()));
    vector<Cell> cur(vector<Cell>(model.states(), Cell()));
    
    // Initialize first round in prev
    for (size_t i = 0; i < model.states(); i++) {
        prev[i].prob = model.pi[i] + model.emissionProb(observation[0], i);
        prev[i].list = shared_ptr<List>(nullptr);
    }
    
    // Recursion
    for (size_t l = 1; l < observation.length(); l++) {
        char symbol = observation[l];
        
        // Cells with positive prob
        /* for (size_t i = 0; i < model.states(); i++) {
         if (prev[i].prob > -numeric_limits<double>::infinity())
         cout << i << " ";
         }
         cout << endl; */
        
        //if (l % 100000 == 0)
        //    cout << (double) l / observation.length() * 100 << "%" << endl;
        
        for (size_t i = 0; i < model.states(); i++) {
            // Fill out a cell
            
            // Find where we should come from
            pair<int, double> best = make_pair(-1, -numeric_limits<double>::infinity());
            for (auto k : model.incommingStates(i)) {
                pair<int, double> candidate;
                candidate = make_pair(k, prev[k].prob + model.transitionProb(k, i));

                if (candidate.second > best.second)
                    best = candidate;
            }
            
            if (best.first == -1) {
                // State is not possible
                cur[i].prob = -numeric_limits<double>::infinity();
                cur[i].list = shared_ptr<List>(nullptr);
            } else {
                // Update current cell with right values
                cur[i].prob = model.emissionProb(symbol, i) + best.second;
                cur[i].list = shared_ptr<List>(new List(best.first, prev[best.first].list));
            }
        }
        
        swap(prev, cur);
    }
    
    // Final result is now in prev
    pair<size_t, double> best = make_pair(-1, -numeric_limits<double>::infinity());
    for (size_t i = 0; i < model.states(); i++) {
        if (prev[i].prob > best.second)
            best = make_pair(i, prev[i].prob);
    }
    
    // Backtrack
    vector<size_t> stateTrace;
    stateTrace.reserve(observation.length());
    stateTrace.push_back(best.first);
    auto node = prev[best.first].list;
    while (node != nullptr) {
        stateTrace.push_back(node->state);
        node = node->prev;
    }
    
    return make_pair(best.second,
                     vector<size_t>(stateTrace.rbegin(), stateTrace.rend()));
}
