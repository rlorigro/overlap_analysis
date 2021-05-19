#include "OverlapMap.hpp"
#include "Graph.hpp"
#include "graph_utils.hpp"

using overlap_analysis::DoubleStrandedGraph;
using overlap_analysis::load_paf_as_graph;
using overlap_analysis::load_adjacency_csv_as_graph;
using overlap_analysis::RegionalOverlapMap;
using overlap_analysis::EdgeDiff;

#include <iostream>
#include <string>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::absolute;
using std::experimental::filesystem::exists;
using std::experimental::filesystem::path;
using std::runtime_error;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::string;
using std::cerr;
using std::cout;


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


void construct_graph(path file_path, DoubleStrandedGraph& graph, uint32_t min_quality){
    RegionalOverlapMap _;

    if (file_path.extension().string() == ".csv"){
        cerr << "Constructing graph from csv file...\n";
        load_adjacency_csv_as_graph(file_path, graph);
    }
    else if (file_path.extension().string() == ".paf"){
        cerr << "Constructing graph from paf file...\n";
        load_paf_as_graph(file_path, _, graph, min_quality);
    }
    else{
        throw runtime_error("ERROR: overlap format does not match 'csv' or 'paf' extension type");
    }
}


void evaluate_overlaps(
        path path_a,
        path path_b,
        path output_directory,
        uint32_t min_quality,
        path excluded_reads_path,
        uint16_t label_type,
        bool plot,
        bool disable_node_union
        ){

    if (exists(output_directory)){
        throw runtime_error("ERROR: output directory already exists");
    }
    else{
        create_directories(output_directory);
    }

    DoubleStrandedGraph graph_a;
    DoubleStrandedGraph graph_b;

    construct_graph(path_a, graph_a, min_quality);
    construct_graph(path_b, graph_b, min_quality);

    // By default, remove all nodes that arent shared by both graphs
    if (not disable_node_union){
        graph_a.node_union(graph_b);
        graph_b.node_union(graph_a);
    }

    if (not excluded_reads_path.empty()){
        exclude_reads_from_graph(graph_a, excluded_reads_path, output_directory/"a_excluded_reads.txt");
        exclude_reads_from_graph(graph_b, excluded_reads_path, output_directory/"b_excluded_reads.txt");
    }

    if (plot) {
        render_graph(graph_a, label_type, output_directory, "a");
        render_graph(graph_b, label_type, output_directory, "b");
    }

    cerr << "Evaluating edge differences..." << '\n';
    EdgeDiff diff(graph_a, graph_b, output_directory);

    ofstream out_file(output_directory / "diff.txt");
    out_file << diff;
}


int main(int argc, char* argv[]){
    path path_a;
    path path_b;
    path output_directory;
    path excluded_reads_path;
    uint32_t min_quality;
    uint16_t label_type;
    string subgraph_argument;
    string subgraph_node_name;
    uint32_t subgraph_radius;
    bool plot;
    bool disable_node_union;

    options_description options("Arguments:");

    options.add_options()
            ("a",
             value<path>(&path_a)
             ->required(),
             "File path of file from which overlap can be inferred\n:"
             "\tPAF file containing alignments to some reference\n "
             "\tCSV file containing a list of reads pairs, with an indication for whether the overlap is cross-strand\n")

            ("b",
             value<path>(&path_b)
             ->required(),
             "File path of file from which overlap can be inferred:\n"
             "\tPAF file containing alignments to some reference \n"
             "\tCSV file containing a list of reads pairs, with an indication for whether the overlap is cross-strand\n")

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

            ("plot",
             bool_switch(&plot)->
             default_value(false),
             "Whether to render the graph as SVG or not. \n"
             "Strongly discouraged for graphs with more than a few thousand nodes")

            ("disable_node_union",
             bool_switch(&disable_node_union)->
             default_value(false),
             "Whether to allow nodes that aren't in both graphs. \n"
             "If this option is specified, all nodes that aren't shared in a and b graph will NOT be deleted")
            ;

    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cerr << options << "\n";
        return 0;
    }
    notify(vm);

    evaluate_overlaps(
            path_a,
            path_b,
            output_directory,
            min_quality,
            excluded_reads_path,
            label_type,
            plot,
            disable_node_union);

    return 0;
}
