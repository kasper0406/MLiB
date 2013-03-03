#include <map>
#include <set>
#include <vector>
#include <string>
#include <memory>

#include "CountingTrainer.h"
#include "Matrix.h"

using namespace std;

unique_ptr<HMM> train_by_counting(vector<string> observations, vector<vector<string>> annotations)
{
    int N = 0, S = 0;
    
    vector<string> states;
    map<string, size_t> translateState;
    vector<char> symbols;
    map<char, size_t> translateSymbol;
    
    for (int i = 0; i < observations.size(); i++) {
        for (string s : annotations[i]) {
            if (translateState.count(s) == 0) {
                translateState.insert(make_pair(s, N++));
                states.push_back(s);
            }
        }
        
        for (char c : observations[i]) {
            if (translateSymbol.count(c) == 0) {
                translateSymbol.insert(make_pair(c, S++));
                symbols.push_back(c);
            }
        }
    }
    
    map<size_t, map<char, size_t>> emissions;
    Matrix<size_t> A(N, N, 0);
    Matrix<size_t> phi(N, S, 0);
    vector<size_t> pi(N, 0);
    
    for (int i = 0; i < observations.size(); i++) {
        for (int j = 0; j < observations[i].length()-1; j++) {
            // Count transition
            A(translateState.at(annotations[i][j]), translateState.at(annotations[i][j+1]))++;
            
            // Count emission
            phi(translateState.at(annotations[i][j]), translateSymbol.at(observations[i][j]))++;
        }
        phi(translateState.at(annotations[i][observations[i].length()-1]),
            translateSymbol.at(observations[i][observations[i].length()-1]))++;
        pi[translateState.at(annotations[i][0])]++;
    }
    
    auto sumA = [N, &A] (size_t row) {
        size_t sum = 0;
        for (size_t j = 0; j < N; j++)
            sum += A(row, j);
        return sum;
    };
    
    auto sumPhi = [S, &phi] (size_t row) {
        size_t sum = 0;
        for (size_t j = 0; j < S; j++)
            sum += phi(row, j);
        return sum;
    };
    
    auto sumPi = [N, &pi] () {
        size_t sum = 0;
        for (size_t i = 0; i < N; i++)
            sum += pi[i];
        return sum;
    };
    
    unique_ptr<HMM> model(new HMM(states, symbols));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            model->A(i, j) = (double) A(i, j) / sumA(i);
        
        for (int j = 0; j < S; j++)
            model->phi(i, j) = (double) phi(i, j) / sumPhi(i);
        
        model->pi[i] = (double) pi[i] / sumPi();
    }
    model->logTransform();
    
    return model;
}
