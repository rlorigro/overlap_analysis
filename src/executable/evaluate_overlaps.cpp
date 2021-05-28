#include "OverlapMap.hpp"
#include "Graph.hpp"
#include "graph_utils.hpp"

using overlap_analysis::DoubleStrandedGraph;
using overlap_analysis::AdjacencyMap;
using overlap_analysis::add_paf_edges_to_adjacency_map;
using overlap_analysis::load_adjacency_csv_as_graph;
using overlap_analysis::for_each_edge_in_shasta_adjacency_csv;
using overlap_analysis::write_edges_to_csv;
using overlap_analysis::PafElement;
using overlap_analysis::RegionalOverlapMap;
using overlap_analysis::ShastaLabel;

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


void construct_graph(path file_path, DoubleStrandedGraph& graph, uint32_t min_quality){
    if (file_path.extension().string() == ".csv"){
        cerr << "Constructing graph from csv file...\n";
        load_adjacency_csv_as_graph(file_path, graph);
    }
    else if (file_path.extension().string() == ".paf"){
        cerr << "Constructing graph from paf file...\n";
        load_paf_as_graph(file_path, graph, min_quality);
    }
    else{
        throw runtime_error("ERROR: overlap format does not match 'csv' or 'paf' extension type");
    }
}


void construct_adjacency_map(path file_path, uint32_t min_quality, AdjacencyMap& adjacency){
    if (file_path.extension().string() == ".csv"){
        cerr << "Constructing graph from csv file...\n";
        load_adjacency_csv_as_adjacency_map(file_path, adjacency);
    }
    else if (file_path.extension().string() == ".paf"){
        cerr << "Constructing graph from paf file...\n";
        load_paf_as_adjacency_map(file_path, adjacency, min_quality);
    }
    else{
        throw runtime_error("ERROR: overlap format does not match 'csv' or 'paf' extension type");
    }
}


void add_csv_edges_to_adjacency_map(path adjacency_path, AdjacencyMap& adjacency_map){
    for_each_edge_in_shasta_adjacency_csv(adjacency_path, [&](
            string& name_a,
            string& name_b,
            bool is_cross_strand,
            ShastaLabel& label){

        auto id_a = adjacency_map.insert_node(name_a);
        auto id_b = adjacency_map.insert_node(name_b);

        adjacency_map.edges.resize(adjacency_map.id_vs_name.left.size());

        // Create a label that only indicates ref membership
        ShastaLabel ref_label(false, false, false, true);
        adjacency_map.insert_edge(id_a, id_b, is_cross_strand, ref_label);
    });
}


void add_reference_edges_to_adjacency_map(path file_path, uint32_t min_quality, AdjacencyMap& adjacency){
    if (file_path.extension().string() == ".csv"){
        cerr << "Constructing graph from csv file...\n";
        add_csv_edges_to_adjacency_map(file_path, adjacency);
    }
    else if (file_path.extension().string() == ".paf"){
        cerr << "Constructing graph from paf file...\n";
        add_paf_edges_to_adjacency_map(file_path, min_quality, adjacency);
    }
    else{
        throw runtime_error("ERROR: overlap format does not match 'csv' or 'paf' extension type");
    }
}


void evaluate_overlaps(
        path ref_overlap_path,
        path overlap_path,
        path output_directory,
        uint32_t min_quality,
        path excluded_reads_path,
        uint16_t label_type,
        bool plot,
        bool no_intersection
        ){

    if (exists(output_directory)){
        throw runtime_error("ERROR: output directory already exists");
    }
    else{
        create_directories(output_directory);
    }

    AdjacencyMap adjacency;

    construct_adjacency_map(overlap_path, min_quality, adjacency);

    add_reference_edges_to_adjacency_map(ref_overlap_path, min_quality, adjacency);

    if (not excluded_reads_path.empty()){
        // TODO: add chimera/palindrome labeling
    }

    if (plot) {
        // TODO: add subgraph extraction and plotting
    }

    path edge_csv_file_path = output_directory / "labeled_candidates.csv";
    write_edges_to_csv(edge_csv_file_path, adjacency);
}


int main(int argc, char* argv[]){
    path ref_overlap_path;
    path overlap_path;
    path output_directory;
    path excluded_reads_path;
    uint32_t min_quality;
    uint16_t label_type;
    bool plot;
    bool no_intersection;

    options_description options("Arguments:");

    options.add_options()
            ("a",
             value<path>(&ref_overlap_path)
             ->required(),
             "REFERENCE GRAPH: path of file from which overlap can be inferred:\n"
             "\tPAF file containing alignments to some reference\n "
             "\tCSV file containing a list of reads pairs, with an indication for whether the overlap is cross-strand\n")

            ("alignment_directory",
             value<path>(&overlap_path)
             ->required(),
             "Alignments TO BE EVALUATED: path of directory containing Shasta alignment details (one file per alignment)\n")

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
             "Simple txt list of read names (one per line) to be excluded from the graph")

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

            ("no_intersection",
             bool_switch(&no_intersection)->
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
            ref_overlap_path,
            overlap_path,
            output_directory,
            min_quality,
            excluded_reads_path,
            label_type,
            plot,
            no_intersection);

    return 0;
}
