#include "OverlapMap.hpp"
#include "Graph.hpp"
#include "graph_utils.hpp"

using overlap_analysis::DoubleStrandedGraph;
using overlap_analysis::load_paf_as_graph;
using overlap_analysis::load_adjacency_csv_as_graph;
using overlap_analysis::RegionalOverlapMap;
using overlap_analysis::GraphDiff;

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
        uint16_t label_type){

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

    render_graph(graph_a, label_type, output_directory, "a");
    render_graph(graph_b, label_type, output_directory, "b");

    GraphDiff diff(graph_a, graph_b);
    cerr << diff << '\n';
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
            label_type);

    return 0;
}
