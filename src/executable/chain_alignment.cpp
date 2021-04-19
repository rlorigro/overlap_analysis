#include "FastqIterator.hpp"
#include "SvgPlot.hpp"
#include "Kernel.hpp"

using overlap_analysis::FastqElement;
using overlap_analysis::FastqIterator;


#include <algorithm>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <iostream>

using std::max_element;
using std::string;
using std::pow;
using std::vector;
using std::cerr;


#include "boost/program_options.hpp"
#include "boost/bimap.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::bool_switch;
using boost::program_options::value;
using boost::geometry::model::point;
using boost::geometry::model::box;
using boost::geometry::cs::cartesian;
using boost::geometry::index::rtree;
using boost::geometry::index::quadratic;
using boost::geometry::index::intersects;


#include "edlib.h"
#include "mummer/sparseSA.hpp"

using mummer::mummer::sparseSA;
using mummer::mummer::match_t;


void plot_mummer_matches(const vector<match_t>& matches, SvgPlot& plot, size_t line_width){
    string type = "rect";
    string color = "#8115A3";

    for (const auto& m: matches){
        for (size_t i=0; i<size_t(m.len); i+=50) {
            plot.add_point(
                    m.ref + i,
                    m.query + i,
                    type,
                    line_width,
                    color);
        }
    }
}


template <class T> void plot_diagonal_scores(vector<T>& diagonal_scores, SvgPlot& plot){
    for (T x=0; x<diagonal_scores.size(); x++){
        plot.add_line(x,T(0),x,diagonal_scores[x],1,"blue");
    }
}


size_t get_diagonal(size_t y_size, size_t x, size_t y){
    return (y_size - 1) - y + x;
}


size_t compute_all_vs_all(vector <FastqElement>& sequences){
    size_t n_comparisons = 0;

    path output_directory = "plots/";

    if (not exists(output_directory)){
        create_directories(absolute(output_directory));
    }

    GaussianSmoother smoother(5,2);

    for (const auto& s: sequences) {
        for (const auto& s2: sequences) {
            if (s.name == s2.name) {
                continue;
            }

            const string& ref = s.sequence;
            const string& query = s2.sequence;

            // Set up plot
            path dotplot_path = output_directory / (s.name + "_vs_" + s2.name + ".svg");
            SvgPlot dotplot(dotplot_path, 800, 800, 0, ref.size(), 0, query.size(), true);

            // Construct suffix array
            auto matcher = sparseSA::create_auto(ref.c_str(), ref.size(), 0, true);
            vector<match_t> matches;
            vector<size_t> diagonal_scores(ref.size() + query.size() - 1, 0);
            vector<size_t> smooth_diagonal_scores(ref.size() + query.size() - 1, 0);

            // Search suffix array for exact matches
            matcher.findMAM_each(query, 10, false, [&](const match_t& m){
                matches.emplace_back(m);

                // Accumulate scores for each diagonal
                auto d = get_diagonal(query.size(), m.ref, m.query);

                diagonal_scores[d] += pow(m.len,2);
            });

            cerr << s.name << ' ' << s2.name << ' ' << matches.size() << '\n';

            plot_mummer_matches(matches, dotplot, ref.size()/1000);

            smoother.smooth(diagonal_scores, smooth_diagonal_scores);

            auto max_score = *max_element(diagonal_scores.begin(), diagonal_scores.end());
            cerr << max_score << '\n';

            path diagonal_plot_path = output_directory / (s.name + "_vs_" + s2.name + "_diagonals.svg");
            SvgPlot diagonal_plot(diagonal_plot_path, 800, 800, 0, ref.size() + query.size(), 0, max_score, true);

            plot_diagonal_scores(smooth_diagonal_scores, diagonal_plot);

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

    while (fastq_iterator.next_fastq_element(element)){
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
             "File path of PAF file containing alignments to some reference")
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

