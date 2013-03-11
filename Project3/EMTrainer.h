#pragma once

#include <vector>
#include <string>

#include "HMM.h"

using namespace std;

void train_by_baumwelch(HMM& model, vector<string> observations);
