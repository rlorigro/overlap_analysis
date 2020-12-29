#ifndef OVERLAP_ANALYSIS_ALIGNMENTINTERVAL_HPP
#define OVERLAP_ANALYSIS_ALIGNMENTINTERVAL_HPP

#include "boost/icl/split_interval_map.hpp"
#include "boost/icl/interval_map.hpp"
#include "boost/icl/interval.hpp"

using boost::icl::partial_enricher;
using boost::icl::total_enricher;
using boost::icl::partial_absorber;
using boost::icl::total_absorber;
using boost::icl::split_interval_map;
using boost::icl::interval_map;
using boost::icl::interval;

#include <unordered_map>
#include <set>
#include <string>

using std::unordered_map;
using std::set;
using std::string;


class RegionalOverlapMap {
public:
    unordered_map <string, interval_map <uint32_t, set<uint32_t>, total_enricher> > intervals;
    uint32_t count;

    void insert(string& region_name, uint32_t start, uint32_t stop);

    RegionalOverlapMap();
};


#endif //OVERLAP_ANALYSIS_ALIGNMENTINTERVAL_HPP
