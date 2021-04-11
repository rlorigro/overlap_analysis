//
//  wavefront.hpp
//  
//
//  Created by Jordan Eizenga on 4/4/21.
//

#ifndef OVERLAP_ANALYSIS_WAVEFRONT_HPP
#define OVERLAP_ANALYSIS_WAVEFRONT_HPP

#include <cstdint>
#include <vector>
#include <string>
#include <limits>

#include <experimental/filesystem>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::absolute;
using std::experimental::filesystem::exists;
using std::experimental::filesystem::path;


/*
 * Score parameters for wavefront alignment. On opening a insertion or
 * deletion, *both* the gap open and gap extend penalties are applied.
 */
struct WFScores {
    WFScores(uint32_t mismatch, uint32_t gap_open, uint32_t gap_extend) :
            mismatch(mismatch), gap_open(gap_open), gap_extend(gap_extend) {}
    uint32_t mismatch;
    uint32_t gap_open;
    uint32_t gap_extend;
};

/*
 * Operation in a CIGAR string.
 */
struct CIGAROp {
    CIGAROp(char op, uint32_t len) : op(op), len(len) {}
    char op;
    uint32_t len;
};

/// Align two sequences using the wavefront alignment algorithm, returns a CIGAR
/// string. Max score is computed as a * (furthest antidiagonal) + b. Returns an empty
/// vector if at any point in the iteration the minimum possible score increases above
/// the max score.
/// Optionally performs pruning (if prune_diff >= 0) of diagonals that are the
/// indicated difference behind the leading diagonal, measured in antidiagonals.
std::vector<CIGAROp> wavefront_align(const std::string& seq1, const std::string& seq2,
                                     const WFScores& scores,
                                     double max_score_a = 0.0,
                                     int32_t max_score_b = std::numeric_limits<int32_t>::max(),
                                     int32_t prune_diff = -1,
                                     path viz_path = "");



#endif /* OVERLAP_ANALYSIS_WAVEFRONT_HPP */
