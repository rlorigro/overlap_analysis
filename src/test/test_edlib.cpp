#include "FastqIterator.hpp"
#include "SvgPlot.hpp"

using overlap_analysis::FastqElement;
using overlap_analysis::FastqIterator;

#include <chrono>
#include <string>
#include <vector>
#include <iostream>

using std::string;
using std::vector;
using std::cerr;

#include "boost/program_options.hpp"
#include "boost/bimap.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::bool_switch;
using boost::program_options::value;

#include "edlib.h"
#include "mummer/sparseSA.hpp"


void plot_edlib_alignment(const string& ref, const string& query, SvgPlot& plot){
    auto config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, nullptr, 0);

    EdlibAlignResult result = edlibAlign(ref.c_str(), ref.size(), query.c_str(), query.size(), config);
    if (result.status != EDLIB_STATUS_OK) {
        throw runtime_error("ERROR: edlib alignment failed");
    }

    vector<size_t> ref_move = {1,1,0,1};
    vector<size_t> query_move = {1,0,1,1};

    size_t prev_ref_index = 0;
    size_t prev_query_index = 0;
    size_t ref_index = ref_move[result.alignment[0]];
    size_t query_index = query_move[result.alignment[0]];

    for (size_t i=0; i<200 and i<size_t(result.alignmentLength); i++) {
        cerr << int(result.alignment[i]);
    }
    cerr << '\n';

    cerr << ref.substr(0,200) << '\n';
    cerr << query.substr(0,200) << '\n';

    for (size_t i=1; i<size_t(result.alignmentLength); i++){
        if (i > 0){
            if (result.alignment[i-1] != result.alignment[i]) {
                plot.add_line(prev_ref_index, prev_query_index, ref_index, query_index, ref.size()/2000, "orange");

                prev_query_index = query_index;
                prev_ref_index = ref_index;
            }
        }

        query_index += query_move[result.alignment[i]];
        ref_index += ref_move[result.alignment[i]];
    }

    plot.add_line(prev_ref_index, prev_query_index, ref_index, query_index, ref.size()/2000, "orange");

    edlibFreeAlignResult(result);
}


void plot_mummer_matches(const string& ref, const string& query, SvgPlot& plot){
    auto matcher = mummer::mummer::sparseSA::create_auto(ref.c_str(), ref.size(), 0, true);

    vector<mummer::mummer::match_t> mams;

    string type = "rect";
    string color = "#8115A3";

    matcher.findMAM_each(query, 12, false, [&](const mummer::mummer::match_t& match){
        mams.emplace_back(match);

        for (long i=0; i<match.len; i+=50) {
            plot.add_point(
                    match.ref + i,
                    match.query + i,
                    type,
                    ref.size()/2000,
                    color);
        }
    });
}


size_t compute_all_vs_all(vector <FastqElement>& sequences){
    size_t n_comparisons = 0;

    path output_directory = "plots/";

    if (not exists(output_directory)){
        create_directories(absolute(output_directory));
    }

    for (const auto& s: sequences) {
        for (const auto& s2: sequences) {
            if (s.name == s2.name) {
                continue;
            }

            const string& ref = s.sequence;
            const string& query = s2.sequence;

            path plot_path = output_directory / (s.name + "_vs_" + s2.name + ".svg");
            SvgPlot plot(plot_path, 800, 800, 0, ref.size(), 0, query.size(), true);

            plot_mummer_matches(ref, query, plot);
            plot_edlib_alignment(ref, query, plot);

            n_comparisons++;
        }
    }

    return n_comparisons;
}


void test(path fastq_path){
    // Pre-load all the sequences into memory using a vector
    vector <FastqElement> sequences;
    FastqElement element;

    FastqIterator fastq_iterator(fastq_path);

    cerr << "Loading sequences..." << '\n';

    while (fastq_iterator.next_element(element)){
        sequences.emplace_back(element);
        cerr << "Read '" << element.name << "' = " << element.sequence.size() << " bp " <<'\n';
    }

    compute_all_vs_all(sequences);
}


int main(int argc, char* argv[]){
    path fastq_path;
    path output_directory;

    options_description options("Arguments:");

    options.add_options()
            ("fastq,i",
             value<path>(&fastq_path)
                     ->required(),
             "File path of Fastq sequences to be aligned all-vs-all")
            ;

    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cerr << options << "\n";
        return 0;
    }
    notify(vm);

    test(fastq_path);

    return 0;
}

