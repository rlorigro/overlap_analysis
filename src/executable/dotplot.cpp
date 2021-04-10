#include "FastqIterator.hpp"
#include "Plot.hpp"
#include "mummer/sparseSA.hpp"

using overlap_analysis::FastqElement;
using overlap_analysis::FastqIterator;

#include <experimental/filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::path;
using std::runtime_error;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::string;
using std::vector;
using std::cerr;
using std::cout;

#include "boost/program_options.hpp"
#include "boost/bimap.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::bool_switch;
using boost::program_options::value;


char complement(char c){
    c = char(toupper(c));

    if (c == 'A'){
        return 'T';
    }
    else if (c == 'C'){
        return 'G';
    }
    else if (c == 'G'){
        return 'C';
    }
    else if (c == 'T'){
        return 'A';
    }
    else {
        throw runtime_error("ERROR: non DNA character cannot be complemented: " + c);
    }
}


void reverse_complement(string& sequence, string& rc_sequence){
    rc_sequence.resize(0);

    for (std::string::reverse_iterator iter=sequence.rbegin(); iter!=sequence.rend(); ++iter){
        rc_sequence += complement(*iter);
    }
}


void dotplot(path fastq_path, path output_directory){
    uint32_t min_length = 12;

    if (exists(output_directory)){
        throw runtime_error("ERROR: output directory already exists");
    }
    else{
        create_directories(absolute(output_directory));
    }

    FastqElement element;
    string rc_sequence;

    FastqIterator fastq_iterator(fastq_path);

    while (fastq_iterator.next_fastq_element(element)){
        path plot_path = output_directory / (element.name + "_dotplot.svg");
        Plot plot(plot_path, 800, 800, 18, 2);

        vector <array <int64_t,2> > coords;

        reverse_complement(element.sequence, rc_sequence);

        if (element.sequence.size() < 400) {
            cerr << element.sequence << '\n';
            cerr << rc_sequence << '\n';
        }

        auto matcher = mummer::mummer::sparseSA::create_auto(element.sequence.c_str(), element.sequence.size(), 0, true);

        vector<mummer::mummer::match_t> mams;

        // min_len must be > 1
        matcher.findMAM_each(rc_sequence, min_length, false, [&](const mummer::mummer::match_t& match){
            mams.emplace_back(match);

            for (int64_t i=0; i<match.len; i+=50){
                array <int64_t,2> coord = {match.ref + i, match.query + i};

                coords.emplace_back(coord);
            }
        });

        plot.add_points(coords, "", 0.1, "");
        plot.generate();
    }
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

            ("output_directory,o",
             value<path>(&output_directory)
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

    dotplot(fastq_path, output_directory);

    return 0;
}




