#include "OverlapMap.hpp"
#include "PafElement.hpp"
#include "Graph.hpp"
#include "graph_utils.hpp"

#include "boost/program_options.hpp"
#include "boost/bimap.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::bool_switch;
using boost::program_options::value;
using boost::bimap;

#include <experimental/filesystem>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <queue>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::absolute;
using std::experimental::filesystem::exists;
using std::experimental::filesystem::path;
using std::runtime_error;
using std::unordered_map;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::string;
using std::vector;
using std::queue;
using std::pair;
using std::cerr;
using std::cout;

typedef bimap<uint32_t,string> uint32_string_bimap;
typedef uint32_string_bimap::value_type bimap_pair;


namespace overlap_analysis {

/// Given an overlap map with the structure: [interval_start, interval_stop) -> {read_id_0, read_id_1, ... },
/// build the edges of a graph, one edge for each inferred overlap.
/// Graph must have existing nodes, stored by the read ID in a vector.
///
/// Method:
/// Iterate the sets and add any edges that don't exist yet in the graph
///
/// Case 1:
/// s1 = {a,b,c}
/// s2 = {a,b,c,d}
///
/// s2 - s1 = {d}
/// add all edges from (s2 - s1) -> s2
///
///
/// Case 2:
/// s1 = {a,b,c,d}
/// s2 = {a,b,c}
///
/// s2 - s1 = {}
/// Do nothing
///
///
/// Case 3:
/// s1 = {a,b,c}
/// s2 = {a,b,d}
///
/// s2 - s1 = {d}
/// add all edges from (s2 - s1) -> s2
///
void create_graph_edges_from_overlap_map(
        RegionalOverlapMap& overlap_map,
        DoubleStrandedGraph& graph) {

    set<uint32_t> empty_set = {};

    for (auto& item: overlap_map.intervals){
        const auto& overlaps = item.second;
        auto& prev_read_set = empty_set;

        for (auto i = begin(overlaps), e = end(overlaps); i!=e; ++i){
//            const auto& interval = i->first;
            auto& read_set = i->second;

//            cerr << interval << " -> ";
//            for (const auto& id: read_set){
//                cerr << id << " " ;
//            }
//            cerr << '\n';

            for (const auto& id: read_set){
                // If this read id is not in the previous set, it indicates that more edges need to be built
                if (prev_read_set.count(id) == 0){
                    for (const auto& other_id: read_set){
                        if (other_id != id){
                            // Will create the edge if it doesn't exist already
                            graph.add_edge(id, other_id);
                        }
                    }
                }
            }
//            cerr << '\n';

            prev_read_set = read_set;
        }
    }
//    cerr << '\n';

}


/// Parse a PAF file: https://github.com/lh3/miniasm/blob/master/PAF.md using the following data:
///
/// Col (1-based)       Data                        Type
/// -------------       ----                        ----
/// 1                   "Query sequence name"       string (assumed to be numeric for this project, using ONT reads)
/// 6                   "Target sequence name"      string
/// 8                   "Target start..."           int
/// 9                   "Target end..."             int
/// 12                  "Mapping quality"           int
///
void load_paf_as_graph(
        path paf_path,
        RegionalOverlapMap& overlap_map,
        DoubleStrandedGraph& graph,
        uint32_t min_quality) {

    ifstream paf_file(paf_path);

    if (not paf_file.good()) {
        throw runtime_error("ERROR: could not open input file: " + paf_path.string());
    }

    string token;
    string region_name;
    string read_name;
    uint32_t start = 0;
    uint32_t stop = 0;
    uint32_t quality = 0;
    bool is_reverse = false;

    uint64_t n_delimiters = 0;
    uint64_t n_lines = 0;
    char c;

    while (paf_file.get(c)){
        if (c == '\t') {
            if (n_delimiters == 0){
                read_name = token;
            }
            else if (n_delimiters == 4){
                is_reverse = (token == "-");
            }
            else if (n_delimiters == 5){
                region_name = token;
            }
            else if (n_delimiters == 7){
                start = stoi(token);
            }
            else if (n_delimiters == 8){
                stop = stoi(token);
            }
            else if (n_delimiters == 11){
                quality = stoi(token);

                if (quality >= min_quality){
                    uint32_t id;
                    uint32_t forward_id;
                    uint32_t reverse_id;

                    id = graph.add_node(read_name);

                    forward_id = 2 * id;
                    reverse_id = 2 * id + 1;

                    if (not is_reverse) {
                        overlap_map.insert(region_name, start, stop, forward_id);
                    } else {
                        overlap_map.insert(region_name, start, stop, reverse_id);
                    }
                }
            }

            token.resize(0);
            n_delimiters++;
        }
        else if (c == '\n'){
            if (n_delimiters < 11){
                throw runtime_error(
                        "ERROR: file provided does not contain sufficient tab delimiters to be PAF on line: " +
                        to_string(n_lines));
            }

            token.resize(0);
            n_delimiters = 0;
            n_lines++;
        }
        else {
            token += c;
        }
    }

    create_graph_edges_from_overlap_map(overlap_map, graph);
}


void load_excluded_read_names_as_set(path excluded_reads_path, set<string>& excluded_reads) {
    ifstream file(excluded_reads_path);
    if (not file.good()){
        throw runtime_error("ERROR: excluded reads file could not be read: " + excluded_reads_path.string());
    }

    string line;

    while (getline(file, line)){
        excluded_reads.emplace(line.substr(0, line.size()));
    }
}


void exclude_reads_from_graph(
        DoubleStrandedGraph& graph,
        path excluded_reads_path,
        path output_path
){
    ofstream file(output_path);
    if (not file.good()){
        throw runtime_error("ERROR: could not write to output: " + output_path.string());
    }

    set<string> excluded_reads;
    load_excluded_read_names_as_set(excluded_reads_path, excluded_reads);

    for (const string& name: excluded_reads){
        graph.remove_node(name);
        file << name << '\n';
    }
}


/// This method is useful for extracting subgraphs from a PAF, is implemented in a slow manner though, because
/// it re-iterates the PAF
/// \param graph
/// \param nodes
/// \param id_vs_name
/// \param double_stranded
/// \param paf_path
void write_all_read_alignments_as_paf(
        Graph& graph,
        path input_paf_path,
        path output_paf_path) {

    set<string> read_names;
    graph.get_all_read_names(read_names);

    string line;
    string read_name;
    ifstream paf_file(input_paf_path);
    ofstream output_file(output_paf_path);
    while (getline(paf_file, line)){
        for (auto& c: line){
            // Parse the line up to the end of the first token in the PAF (TSV)
            if (c == '\t'){

                // If this read is contained in the graph or subgraph, then write the line to the output PAF
                if (read_names.count(read_name) > 0){
                    output_file << line << '\n';
                }

                // Reset token
                read_name.resize(0);

                // Stop parsing line
                break;
            }
            read_name += c;
        }
    }
}


void write_arguments(
        path file_path,
        path paf_path,
        uint32_t min_quality,
        path excluded_reads_path,
        uint16_t label_type,
        string subgraph_node_name,
        uint32_t subgraph_radius){

    ofstream file(file_path);

    file << "paf_path," << paf_path << '\n';
    file << "min_quality," << min_quality << '\n';
    file << "excluded_reads_path," << excluded_reads_path << '\n';
    file << "label_type," << label_type << '\n';
    file << "subgraph_node_name," << subgraph_node_name << '\n';
    file << "subgraph_radius," << subgraph_radius << '\n';
}


void plot_graph(
        path paf_path,
        path output_directory,
        uint32_t min_quality,
        path excluded_reads_path,
        uint16_t label_type,
        string subgraph_node_name,
        uint32_t subgraph_radius) {

    if (not exists(paf_path)){
        throw runtime_error("ERROR: input PAF does not exist: " + paf_path.string());
    }
    if (not exists(excluded_reads_path) and not excluded_reads_path.empty()){
        throw runtime_error("ERROR: input excluded reads file does not exist: " + excluded_reads_path.string());
    }
    if (exists(output_directory)){
        throw runtime_error("ERROR: output directory already exists");
    }
    else{
        create_directories(output_directory);
    }

    path args_path = output_directory / "args.csv";

    write_arguments(
            absolute(args_path),
            absolute(paf_path),
            min_quality,
            absolute(excluded_reads_path),
            label_type,
            subgraph_node_name,
            subgraph_radius);

    RegionalOverlapMap overlap_map;
    DoubleStrandedGraph overlap_graph;

    cerr << "Loading PAF as graph...\n";
    load_paf_as_graph(paf_path, overlap_map, overlap_graph, min_quality);

    // Delete any nodes from the graph if a list of reads to be excluded is provided
    if (not excluded_reads_path.empty()){
        cerr << "Excluding reads...\n";
        path reads_excluded_log_path = output_directory / "reads_excluded.txt";
        exclude_reads_from_graph(overlap_graph, excluded_reads_path, reads_excluded_log_path);
    }

    // Do subgraph extraction if a node was provided as a start point for BFS
    if (not subgraph_node_name.empty()){
        cerr << "Extracting subgraph of radius " << subgraph_radius << " around node " << subgraph_node_name << '\n';

        path subgraph_paf_path = output_directory / "subgraph.paf";

        Graph subgraph;

        overlap_graph.create_subgraph(subgraph_node_name, subgraph_radius, subgraph);

        // If the user specified a subgraph then it may be useful to write out the subgraph PAF as well
        write_all_read_alignments_as_paf(
                subgraph,
                paf_path,
                subgraph_paf_path);

        render_graph(subgraph, label_type, output_directory);
    }
    else {
        render_graph(overlap_graph, label_type, output_directory);
    }
}

}

using overlap_analysis::plot_graph;

pair<string, uint32_t> parse_subgraph_argument(string subgraph_argument) {
    char delimiter = ':';
    string token;
    string name;
    uint32_t radius;

    size_t n_delimiters = 0;

    for (size_t i=0; i<subgraph_argument.size(); i++){
        char c = subgraph_argument[i];

        if (c == delimiter) {
            if (n_delimiters == 0) {
                name = token;
                token.clear();
                n_delimiters++;
            }
            else{
                throw runtime_error("ERROR: subgraph argument has too many ':' delimiters: " + subgraph_argument);
            }
        }
        else {
            token += c;
        }
    }

    if (n_delimiters != 1){
        throw runtime_error("ERROR: subgraph argument has no ':' delimiters: " + subgraph_argument);
    }

    radius = stoi(token);

    return {name, radius};
}


int main(int argc, char* argv[]){
    path paf_path;
    path output_directory;
    path excluded_reads_path;
    uint32_t min_quality;
    uint16_t label_type;
    string subgraph_argument;
    string subgraph_node_name;
    uint32_t subgraph_radius;

    options_description options("Arguments:");

    options.add_options()
            ("paf_path",
             value<path>(&paf_path)
             ->required(),
             "File path of PAF file containing alignments to some reference")

            ("output_dir",
             value<path>(&output_directory)
             ->required(),
             "Where to dump output SVG, CSV, PAF, etc")

            ("mapq",
            value<uint32_t>(&min_quality)
            ->required()
            ->default_value(50),
            "Minimum allowed mapping quality to load into the graph")

            ("exclude",
             value<path>(&excluded_reads_path)
             ->default_value(""),
             "File path of PAF file containing alignments to some reference")

            ("label,l",
             value<uint16_t>(&label_type)
             ->default_value(0),
             "Labelling scheme for nodes:\n"
             "\t0 - None\n"
             "\t1 - Original names\n"
             "\t2 - Numeric")

            ("subgraph",
             value<string>(&subgraph_argument)
             ->default_value(""),
             "Subset the graph by radius r around a node n, using the format n:r")
            ;


    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cerr << options << "\n";
        return 0;
    }
    notify(vm);

    if (not subgraph_argument.empty()) {
        tie(subgraph_node_name, subgraph_radius) = parse_subgraph_argument(subgraph_argument);
    }
    else{
        subgraph_node_name = "";
        subgraph_radius = 0;
    }

    plot_graph(paf_path, output_directory, min_quality, excluded_reads_path, label_type, subgraph_node_name, subgraph_radius);

    return 0;
}
