#include "AlignmentInterval.hpp"
#include "boost/program_options.hpp"
#include "boost/bimap.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using boost::bimap;

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/graphics.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/fileformats/GraphIO.h>

using ogdf::Graph;
using ogdf::GraphAttributes;
using ogdf::node;
using ogdf::edge;
using ogdf::FMMMLayout;
using ogdf::FMMMOptions;
using ogdf::Shape;

#include <experimental/filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::experimental::filesystem::path;
using std::runtime_error;
using std::ifstream;
using std::to_string;
using std::string;
using std::vector;
using std::cerr;
using std::cout;

typedef bimap<uint32_t,uint32_t> uint32_bimap;
typedef uint32_bimap::value_type bimap_pair;


uint32_t load_csv_as_id_map(path& read_csv_path, uint32_bimap& id_vs_name){
    ///
    /// Parse a line of the ReadSummary.csv with the following format:
    ///     Id,Name, ... etc
    ///     0,512, ... etc
    ///
    /// And return the max id observed
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
    uint32_t max_id = 0;
    char c;

    while (read_info_file.get(c)) {
        if (c == ',') {
            if (n_lines == 0){
                continue;
            }

            if (n_delimiters == 0) {
                id = stoi(token);

                if (id > max_id){
                    max_id = id;
                }
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

    return max_id;
}


/// Parse a PAF file: https://github.com/lh3/miniasm/blob/master/PAF.md using the following data:
///
/// Col (1-based)       Data                        Type
/// 1                   "Query sequence name"       string (assumed to be numeric for this project, using ONT reads)
/// 6                   "Target sequence name"      string
/// 8                   "Target start..."           int
/// 9                   "Target end..."             int
///
void load_paf_as_graph(
        path paf_path,
        uint32_bimap& id_vs_name,
        RegionalOverlapMap& overlap_map,
        ogdf::Graph& graph,
        vector<node>& nodes){

    ifstream paf_file(paf_path);

    if (not paf_file.good()){
        throw runtime_error("ERROR: could not open input file: " + paf_path.string());
    }

    string token;
    string region_name;
    uint32_t read_name;     // Assume numeric read names
    uint32_t start;
    uint32_t stop;
    uint64_t n_delimiters = 0;
    uint64_t n_lines = 0;
    char c;

    while (paf_file.get(c)) {
        if (c == '\t') {
            if (n_delimiters == 0) {
                read_name = stoi(token);
            }
            else if (n_delimiters == 5) {
                region_name = token;
            }
            else if (n_delimiters == 7) {
                start = stoi(token);
            }
            else if (n_delimiters == 8) {
                stop = stoi(token);
                cerr << region_name << " " << start << " " << stop << '\n';

                uint32_t id = id_vs_name.right.at(read_name);
                overlap_map.insert(region_name, start, stop, id);
                nodes[id] = graph.newNode();
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

    // TODO: Move this into separate function and rename this existing one as load_as_interval_map or w.e.
    for (auto& item: overlap_map.intervals){
        const auto& overlaps = item.second;
        const auto& prev = item.second.begin();

        // TODO: Build first set of edges by iterating all combos
//        for (){
//
//        }

        // TODO: Iterate the sets and add any edges that don't exist yet in the graph, as they are encountered.
        // Case 1:
        // s1 = {a,b,c}
        // s2 = {a,b,c,d}
        //
        // s2 - s1 = {d}
        // s1 AND s2 = {a,b,c}
        // add all edges from (s2 - s1) -> (s2 AND s1)
        //
        //
        // Case 2:
        // s1 = {a,b,c,d}
        // s2 = {a,b,c}
        //
        // s2 - s1 = {}
        // s1 AND s2 = {a,b,c}
        // Do nothing
        //
        //
        // Case 3:
        // s1 = {a,b,c}
        // s2 = {a,b,d}
        //
        // s2 - s1 = {d}
        // s1 AND s2 = {a,b}
        // add all edges from (s2 - s1) -> (s2 AND s1)
        //

        for (auto i = ++begin(overlaps), e = end(overlaps); i !=e; ++i){
            const auto& interval = i->first;
            const auto& read_set = i->second;
            cerr << interval << " -> ";
            for (auto& item2: read_set){
                cerr << item2 << " " ;
            }
            cerr << '\n';
        }
    }
    cerr << '\n';
}


void label_alignment_info(path info_csv_path, path read_csv_path, path paf_path){
    uint32_bimap id_vs_name;
    uint32_t max_id;
    max_id = load_csv_as_id_map(read_csv_path, id_vs_name);

    RegionalOverlapMap overlap_map;
    Graph overlap_graph;
    vector<node> nodes;
    nodes.resize(max_id+1, nullptr);
    load_paf_as_graph(paf_path, id_vs_name, overlap_map, overlap_graph, nodes);

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
