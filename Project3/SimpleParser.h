#pragma once

#include <vector>
#include <string>

using namespace std;

class SimpleParser
{
public:
    string handleNoncoding(char obs) { return "N"; }
    
    vector<string> handleStart(string obs) { return { "S" }; };
    vector<string> handleEnd(string obs) { return { "E" }; };
    vector<string> handleCoding(string obs) { return { "C" }; };
    
    vector<string> handleReverseStart(string obs) { return { "RS" }; };
    vector<string> handleReverseEnd(string obs) { return { "RE" }; };
    vector<string> handleReverseCoding(string obs) { return { "RC" }; };
};
