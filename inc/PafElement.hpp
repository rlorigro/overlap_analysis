#ifndef OVERLAP_ANALYSIS_PAFELEMENT_HPP
#define OVERLAP_ANALYSIS_PAFELEMENT_HPP

#include <string>

using std::string;

namespace overlap_analysis {

class PafElement {
public:
    string target_name;
    uint32_t start;
    uint32_t stop;
    uint32_t quality;

    PafElement()=default;
    PafElement(string target_name, uint32_t start, uint32_t stop, uint32_t quality);
};

}

#endif //OVERLAP_ANALYSIS_PAFELEMENT_HPP
