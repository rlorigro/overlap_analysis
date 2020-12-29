#include "AlignmentInterval.hpp"
#include "boost/program_options.hpp"
#include "boost/bimap.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using boost::bimap;

#include <experimental/filesystem>
#include <iostream>
#include <fstream>
#include <string>

using std::experimental::filesystem::path;
using std::runtime_error;
using std::ifstream;
using std::to_string;
using std::string;
using std::cerr;
using std::cout;


typedef bimap<uint32_t,uint32_t> uint32_bimap;
typedef uint32_bimap::value_type bimap_pair;


void load_csv_as_id_map(path& read_csv_path, uint32_bimap& id_vs_name){
    ///
    /// Parse a line of the ReadSummary.csv with the following format:
    ///     Id,Name, ... etc
    ///     0,512, ... etc
    ///
    ifstream read_info_file(read_csv_path);

    if (not read_info_file.good()){
        throw runtime_error("ERROR: could not open input file: " + read_csv_path.string());
    }

    string token;
    uint32_t id;
    uint32_t name;
    uint64_t n_delimiters = 0;
    uint64_t n_lines = 0;
    char c;

    while (read_info_file.get(c)) {
        if (c == ',') {
            if (n_lines == 0){
                continue;
            }

            if (n_delimiters == 0) {
                id = stoi(token);
            }
            else if (n_delimiters == 1) {
                name = stoi(token);
                id_vs_name.insert(bimap_pair(id, name));
            }

            token.resize(0);
            n_delimiters++;
        }
        else if (c == '\n'){
            token.resize(0);
            n_delimiters = 0;
            n_lines++;
        }
        else {
            token += c;
        }
    }
}


void load_paf_as_graph(path paf_path, RegionalOverlapMap& overlap_map){
    /// Parse a PAF file: https://github.com/lh3/miniasm/blob/master/PAF.md
    /// Using the "Target sequence name", "Target start...", and "Target end..." columns (6,8,9)
    ifstream paf_file(paf_path);

    if (not paf_file.good()){
        throw runtime_error("ERROR: could not open input file: " + paf_path.string());
    }

    string token;
    string region_name;
    uint32_t start;
    uint32_t stop;
    uint64_t n_delimiters = 0;
    uint64_t n_lines = 0;
    char c;

    while (paf_file.get(c)) {
        if (c == '\t') {
            if (n_delimiters == 5) {
                region_name = token;
            }
            else if (n_delimiters == 7) {
                start = stoi(token);
            }
            else if (n_delimiters == 8) {
                stop = stoi(token);
                cerr << region_name << " " << start << " " << stop << '\n';
                overlap_map.insert(region_name, start, stop);
            }

            token.resize(0);
            n_delimiters++;
        }
        else if (c == '\n'){
            token.resize(0);
            n_delimiters = 0;
            n_lines++;
        }
        else {
            token += c;
        }
    }
}


void label_alignment_info(path info_csv_path, path read_csv_path, path paf_path){
    uint32_bimap id_vs_name;
    load_csv_as_id_map(read_csv_path, id_vs_name);

    RegionalOverlapMap overlap_map;
    load_paf_as_graph(paf_path, overlap_map);

    for (auto& item: overlap_map.intervals){
        cerr << item.first << '\n';
        for (auto& item1: item.second){
            cerr << item1.first << " -> ";
            for (auto& item2: item1.second){
                cerr << item2 << " " ;
            }
            cerr << '\n';
        }
    }
    cerr << '\n';

    {
        uint32_t i = 27;
        cerr << "finding: " << i << '\n';
        auto iter = overlap_map.intervals.at("A").find(i);
        cerr << iter->first << " -> ";
        for (auto& item: iter->second){
            cerr << item << " " ;
        }
        cerr << '\n';
    }
    {
        uint32_t i = 50;
        cerr << "finding: " << i << '\n';
        auto iter = overlap_map.intervals.at("A").find(i);
        cerr << iter->first << " -> ";
        for (auto& item: iter->second){
            cerr << item << " " ;
        }
        cerr << '\n';
    }
    {
        uint32_t i = 75;
        cerr << "finding: " << i << '\n';
        auto iter = overlap_map.intervals.at("A").find(i);
        cerr << iter->first << " -> ";
        for (auto& item: iter->second){
            cerr << item << " " ;
        }
        cerr << '\n';
    }
}


int main(int argc, char* argv[]){
    path info_csv_path;
    path read_csv_path;
    path paf_path;

    options_description options("Arguments:");

    options.add_options()
            ("info_csv_path",
             value<path>(&info_csv_path),
             "File path of CSV alignment info dump created by shasta readGraph.creationMethod 2 with debug mode on")

            ("read_csv_path",
             value<path>(&read_csv_path),
             "File path of CSV file ReadSummary.csv produced by shasta")

            ("paf_path",
             value<path>(&paf_path),
             "File path of PAF file containing alignments to some reference")
            ;

    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cerr << options << "\n";
        return 0;
    }

    label_alignment_info(info_csv_path, read_csv_path, paf_path);

    return 0;
}
