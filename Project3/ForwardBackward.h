#pragma once

#include <string>
#include <vector>
#include <tuple>

#include "Matrix.h"
#include "HMM.h"

using namespace std;

tuple<vector<double>,Matrix<double>,Matrix<double>> forward_backward(string obs, const HMM& model);
