#include "PafElement.hpp"

namespace overlap_analysis{

PafElement::PafElement(string target_name, uint32_t start, uint32_t stop, uint32_t quality):
    target_name(target_name),
    start(start),
    stop(stop),
    quality(quality)
{}

}