#include "OverlapMap.hpp"
#include "Graph.hpp"
#include "graph_utils.hpp"

using overlap_analysis::DoubleStrandedGraph;
using overlap_analysis::load_paf_as_graph;
using overlap_analysis::load_adjacency_csv_as_graph;
using overlap_analysis::RegionalOverlapMap;
using overlap_analysis::ShastaLabel;
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
        path ref_overlap_path,
        path overlap_path,
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

    DoubleStrandedGraph ref_graph;
    DoubleStrandedGraph graph;

    construct_graph(ref_overlap_path, ref_graph, min_quality);
    construct_graph(overlap_path, graph, min_quality);

    // By default, remove all nodes that arent shared by both graphs
    if (not disable_node_union){
        ref_graph.node_union(graph);
        graph.node_union(ref_graph);
    }

    if (not excluded_reads_path.empty()){
        exclude_reads_from_graph(ref_graph, excluded_reads_path, output_directory / "excluded_ref_overlaps.txt");
        exclude_reads_from_graph(graph, excluded_reads_path, output_directory / "excluded_overlaps.txt");
    }

    if (plot) {
        render_graph(ref_graph, label_type, output_directory, "ref_overlap_graph");
        render_graph(graph, label_type, output_directory, "overlap_graph");
    }

    cerr << "Evaluating edge differences..." << '\n';
    EdgeDiff diff;

    path edge_table_filename = output_directory / "labeled_overlaps.csv";
    ofstream edge_table_file(edge_table_filename);
    if (not edge_table_file.is_open() or not edge_table_file.good()){
        throw runtime_error("ERROR: couldn't write file: " + edge_table_filename.string());
    }

    diff.for_each_edge_comparison(
            ref_graph,
            graph,
            [&](uint32_t id0,
                uint32_t id1,
                bool is_cross_strand,
                bool in_ref,
                bool in_non_ref){

        ShastaLabel label;
        auto success = graph.find_label(id0, id1, is_cross_strand, label);

        // TODO add names to the parameter list and consider using a class to transfer data...
        // id is insufficient for writing the output table
        if (success){
            edge_table_file << "stuff" << '\n';
        }
    });

    ofstream summary_file(output_directory / "summary.txt");
    if (not summary_file.is_open() or not summary_file.good()){
        throw runtime_error("ERROR: couldn't write file: " + edge_table_filename.string());
    }
    summary_file << diff;
}


int main(int argc, char* argv[]){
    path ref_overlap_path;
    path overlap_path;
    path output_directory;
    path excluded_reads_path;
    uint32_t min_quality;
    uint16_t label_type;
    bool plot;
    bool disable_node_union;

    options_description options("Arguments:");

    options.add_options()
            ("a",
             value<path>(&ref_overlap_path)
             ->required(),
             "File path of file from which overlap can be inferred\n:"
             "\tPAF file containing alignments to some reference\n "
             "\tCSV file containing a list of reads pairs, with an indication for whether the overlap is cross-strand\n")

            ("b",
             value<path>(&overlap_path)
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
            ref_overlap_path,
            overlap_path,
            output_directory,
            min_quality,
            excluded_reads_path,
            label_type,
            plot,
            disable_node_union);

    return 0;
}
