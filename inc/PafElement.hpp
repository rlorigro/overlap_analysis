#ifndef OVERLAP_ANALYSIS_PAFELEMENT_HPP
#define OVERLAP_ANALYSIS_PAFELEMENT_HPP

#include <string>
#include <iostream>

using std::string;
using std::ostream;
using std::cerr;

namespace overlap_analysis {


class PafElement {
public:
    string target_name;
    string query_name;
    uint32_t start;
    uint32_t stop;
    uint32_t map_quality;
    bool is_reverse;

    PafElement()=default;
    PafElement(string& target_name,
               string& query_name,
               uint32_t start,
               uint32_t stop,
               uint32_t map_quality,
               bool is_reverse);
};

ostream& operator<<(ostream& o, PafElement& e);

}

#endif //OVERLAP_ANALYSIS_PAFELEMENT_HPP
