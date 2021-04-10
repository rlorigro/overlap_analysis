#ifndef OVERLAP_ANALYSIS_DOTPLOT_HPP
#define OVERLAP_ANALYSIS_DOTPLOT_HPP

#include <experimental/filesystem>
#include <string>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::path;
using std::runtime_error;
using std::string;


void generate_dotplot(
        string& name,
        string& ref_sequence,
        string& query_sequence,
        uint32_t min_length,
        path& output_directory
        );

#endif //OVERLAP_ANALYSIS_DOTPLOT_HPP
