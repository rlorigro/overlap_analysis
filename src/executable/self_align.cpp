#include "FastqIterator.hpp"
#include "wavefront.hpp"

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


void self_align(path fastq_path){
    uint32_t mismatch_penalty = 1;

    WFScores scores(mismatch_penalty, 1, 1);

    double error_rate = 0.15;

    double a = 2 * error_rate * mismatch_penalty;
    int32_t prune = -2 * log(1e-6) / log(4);

    FastqElement element;
    string rc_sequence;

    FastqIterator fastq_iterator(fastq_path);

    while (fastq_iterator.next_fastq_element(element)){
        path visualization_path = element.name + "_matrix.svg";

        int32_t b = 1 * sqrt(element.sequence.size() + element.sequence.size()) * mismatch_penalty;

        cerr << "self-aligning sequence " << element.name << " of length " << element.sequence.size() << " with params: \n"
             << '\t' << "a=" << a << '\n'
             << '\t' << "b=" << b << '\n'
             << '\t' << "prune=" << prune << '\n';

        reverse_complement(element.sequence, rc_sequence);

        if (element.sequence.size() < 400) {
            cerr << element.sequence << '\n';
            cerr << rc_sequence << '\n';
        }

        auto cigars = wavefront_align(element.sequence, rc_sequence, scores, a, b, prune, visualization_path);

        if (cigars.empty()){
            cerr << "alignment aborted" << '\n';
        }
        else {
            for (const auto& item: cigars) {
                cerr << item.len << item.op << ' ';
            }
        }
        cerr << '\n';
    }
}


int main(int argc, char* argv[]){
    path fastq_path;

    options_description options("Arguments:");

    options.add_options()
            ("fastq",
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

    self_align(fastq_path);

    return 0;
}

