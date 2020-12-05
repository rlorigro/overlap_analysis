#ifndef OVERLAP_ANALYSIS_SEQUENCEELEMENT_HPP
#define OVERLAP_ANALYSIS_SEQUENCEELEMENT_HPP

#include <vector>
#include <string>

using std::vector;
using std::string;


class SequenceElement{
public:
    string sequence;
    string name;

    SequenceElement()=default;
    SequenceElement(string name, string sequence);
};

#endif //OVERLAP_ANALYSIS_SEQUENCEELEMENT_HPP
