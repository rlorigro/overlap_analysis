#ifndef OVERLAP_ANALYSIS_DOTPLOT_HPP
#define OVERLAP_ANALYSIS_DOTPLOT_HPP

#include "SvgPlot.hpp"

#include <experimental/filesystem>
#include <string>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::path;
using std::runtime_error;
using std::string;


void reverse_complement(string& sequence, string& rc_sequence);

void generate_dotplot(
        const string& name,
        const string& ref_sequence,
        const string& query_sequence,
        uint32_t min_length,
        SvgPlot& plot,
        bool both_strands = false,
        int64_t x_offset = 0,
        int64_t y_offset = 0
        );

#endif //OVERLAP_ANALYSIS_DOTPLOT_HPP
