#pragma once

#include <memory>
#include <vector>
#include <string>

#include "HMM.h"

using namespace std;

unique_ptr<HMM> train_by_counting(vector<string> observations, vector<vector<string>> annotations);
