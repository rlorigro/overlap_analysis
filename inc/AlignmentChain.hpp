#ifndef OVERLAP_ANALYSIS_ALIGNMENTCHAIN_HPP
#define OVERLAP_ANALYSIS_ALIGNMENTCHAIN_HPP

#include <experimental/filesystem>
#include <ostream>
#include <vector>
#include <string>
#include <utility>
#include <set>
#include <list>
#include <map>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::path;
using std::ostream;
using std::vector;
using std::string;
using std::pair;
using std::set;
using std::list;
using std::map;


class ChainElement {
public:
    string paf_line;
    string ref_name;
    uint32_t ref_start;
    uint32_t ref_stop;
    uint32_t query_start;
    uint32_t query_stop;
    uint32_t ref_length;
    bool is_reverse;

    /// Methods ///
    ChainElement(
            string& paf_line,
            string& ref_name,
            uint32_t ref_start,
            uint32_t ref_stop,
            uint32_t query_start,
            uint32_t query_stop,
            uint32_t ref_length,
            bool is_reverse);

    uint32_t distance_to_end_of_contig() const;

    ChainElement()=default;
};


ostream& operator<<(ostream& o, const ChainElement& e);


class AlignmentChain {
public:
    const string query_name;
    vector <ChainElement> chain;
    const uint32_t max_gap = 50000;
    const uint32_t gap_penalty = 10000;

    /// Methods ///
    AlignmentChain()=default;
    void add(ChainElement& e);
    void sort_chain();
    void split(set <pair <size_t, size_t> >& subchain_bounds, pair <size_t, size_t> bounds = {0,0});
    uint32_t compute_distance(ChainElement& a, ChainElement& b);
};


class AlignmentChains {
public:
    map <string, AlignmentChain> chains;

    /// Methods ///
    AlignmentChains()=default;
    void add_alignment(string line);
    void load_from_paf(path paf_path);
    void split_all_chains();
};


#endif //OVERLAP_ANALYSIS_ALIGNMENTCHAIN_HPP
