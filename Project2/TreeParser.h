#pragma once

#include <string>
#include <tuple>

using namespace std;

class TreeParser {
public:
    static string handleNoncoding(char obs);
    
    static tuple<string,string,string> handleStart(string obs);
    static tuple<string,string,string> handleEnd(string obs);
    static tuple<string,string,string> handleCoding(string obs);
    
    static tuple<string,string,string> handleReverseStart(string obs);
    static tuple<string,string,string> handleReverseEnd(string obs);
    static tuple<string,string,string> handleReverseCoding(string obs);
};
