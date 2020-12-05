#include "RunlengthSequenceElement.hpp"

RunlengthSequenceElement::RunlengthSequenceElement(string name, string sequence, vector<uint16_t> lengths):
    lengths(lengths),
    sequence(sequence),
    name(name)
{}

