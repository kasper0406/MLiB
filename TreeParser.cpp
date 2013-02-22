#include <tuple>
#include <string>

#include "TreeParser.h"

string TreeParser::handleNoncoding(char obs)
{
    return "NC";
}

tuple<string,string,string> TreeParser::handleStart(string obs)
{
    return make_tuple(string("S1").append(obs.substr(0,1)),
                      string("S2").append(obs.substr(1,1)),
                      string("S3").append(obs.substr(2,1)));
}

tuple<string,string,string> TreeParser::handleEnd(string obs)
{
    return make_tuple(string("E1").append(obs.substr(0,1)),
                      string("E2").append(obs.substr(1,1)),
                      string("E3").append(obs.substr(2,1)));
}

tuple<string,string,string> TreeParser::handleCoding(string obs)
{
    return make_tuple(string("T1").append(obs.substr(0,1)),
                      string("T2").append(obs.substr(0,2)),
                      string("T3").append(obs.substr(0,3)));
}

tuple<string,string,string> TreeParser::handleReverseStart(string obs)
{
    return make_tuple(string("RS1").append(obs.substr(0,1)),
                      string("RS2").append(obs.substr(1,1)),
                      string("RS3").append(obs.substr(2,1)));
}

tuple<string,string,string> TreeParser::handleReverseEnd(string obs)
{
    return make_tuple(string("RE1").append(obs.substr(0,1)),
                      string("RE2").append(obs.substr(1,1)),
                      string("RE3").append(obs.substr(2,1)));
}

tuple<string,string,string> TreeParser::handleReverseCoding(string obs)
{
    return make_tuple(string("RT1").append(obs.substr(0,1)),
                      string("RT2").append(obs.substr(0,2)),
                      string("RT3").append(obs.substr(0,3)));
}
