#include <vector>
#include <string>
#include <iostream>

#include "ViterbiTrainer.h"
#include "HMM.h"
#include "Viterbi.h"
#include "CountingTrainer.h"

using namespace std;

void train_by_viterbi(HMM& model, vector<string> observations, unsigned int iterations)
{
    for (unsigned int i = 1; i <= iterations; i++) {
        model.finalize();
        
        vector<vector<string>> annotations;
        for (string observation : observations) {
            vector<size_t> trace;
            tie(ignore, trace) = viterbi(observation, model);
            
            vector<string> annotation;
            annotation.reserve(trace.size());
            for (size_t state : trace)
                annotation.push_back(model.stateLabel(state));
            
            annotations.push_back(annotation);
        }
        
        model.unlock();
        
        train_by_counting(model, observations, annotations);
        
        cout << "Finished iteration #" << i << " in Viterbi training!" << endl;
    }
}
