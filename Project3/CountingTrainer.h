#pragma once

#include <vector>
#include <string>

#include "HMM.h"

using namespace std;

void train_by_counting(HMM& model, vector<string> observations, vector<vector<string>> annotations);
