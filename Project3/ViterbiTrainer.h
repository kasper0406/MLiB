#pragma once

#include <vector>
#include <string>

#include "HMM.h"

using namespace std;

void train_by_viterbi(HMM& model, vector<string> observations, unsigned int iterations);
