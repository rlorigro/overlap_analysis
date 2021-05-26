#include "PafElement.hpp"

namespace overlap_analysis{

PafElement::PafElement(string& target_name,
                       string& query_name,
                       uint32_t start,
                       uint32_t stop,
                       uint32_t map_quality,
                       bool is_reverse):
    target_name(target_name),
    query_name(query_name),
    start(start),
    stop(stop),
    map_quality(map_quality),
    is_reverse(is_reverse)
{}

}