#include "AlignmentInterval.hpp"


// TODO: add some kind of pointer/reference to node in graph (read id?)
void RegionalOverlapMap::insert(string& region_name, uint32_t start, uint32_t stop) {
    // Check if this region has been encountered before, and initialize it if needed
    if (intervals.count(region_name) == 0){
        interval_map <uint32_t, set<uint32_t>, total_enricher> empty_interval_map;
        intervals.emplace(region_name, empty_interval_map);
    }

    // Construct a Boost interval for this numeric range/coord
    auto a = interval<uint32_t>::right_open(start, stop);
    set<uint32_t> b = {++count};

    // Within this contig, add the numeric interval
    intervals.at(region_name).add({a, b});
}


RegionalOverlapMap::RegionalOverlapMap():
    intervals(),
    count(0)
{}
