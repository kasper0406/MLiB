#pragma once

#include <vector>
#include <string>

#include "HMM.h"

using namespace std;

pair<double,vector<size_t>> viterbi(string observation, const HMM& model);
