#ifndef OVERLAP_ANALYSIS_RUNLENGTHSEQUENCEELEMENT_HPP
#define OVERLAP_ANALYSIS_RUNLENGTHSEQUENCEELEMENT_HPP

#include <string>
#include <vector>

using std::vector;
using std::string;


class RunlengthSequenceElement{
public:
    vector<uint16_t> lengths;
    string sequence;
    string name;

    RunlengthSequenceElement()=default;
    RunlengthSequenceElement(string name, string sequence, vector<uint16_t> lengths);
};

#endif //OVERLAP_ANALYSIS_RUNLENGTHSEQUENCEELEMENT_HPP
