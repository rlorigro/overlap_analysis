#include "OverlapMap.hpp"
#include "Graph.hpp"
#include "graph_utils.hpp"

using overlap_analysis::DoubleStrandedGraph;
using overlap_analysis::AdjacencyMap;
using overlap_analysis::add_paf_edges_to_adjacency_map;
using overlap_analysis::load_adjacency_csv_as_graph;
using overlap_analysis::for_each_edge_in_shasta_adjacency_csv;
using overlap_analysis::for_each_edge_in_adjacency;
using overlap_analysis::write_edges_to_csv;
using overlap_analysis::PafElement;
using overlap_analysis::RegionalOverlapMap;
using overlap_analysis::ShastaLabel;

#include <unordered_set>
#include <iostream>
#include <string>

using std::experimental::filesystem::directory_iterator;
using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::absolute;
using std::experimental::filesystem::rename;
using std::experimental::filesystem::exists;
using std::experimental::filesystem::path;
using std::unordered_set;
using std::runtime_error;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::string;
using std::cerr;
using std::cout;


void load_alignment_files_as_adjacency_map(
        path alignment_directory,
        AdjacencyMap& adjacency_map,
        vector <path>& alignment_files){

    for(auto& file_iter: directory_iterator(alignment_directory)){
        string filename = file_iter.path().filename().string();

        if (filename.empty()){
            throw runtime_error("ERROR: filename in alignment directory is empty: " + file_iter.path().string());
        }

        size_t n_delimiters = 0;
        string token;
        string name_a;
        string name_b;
        bool is_same_strand;

        for (char c: filename){
            if (c == '_' or c == '.'){
                if (n_delimiters == 0){
                    name_a = token;
                }
                else if (n_delimiters == 1){
                    name_b = token;
                }
                else if (n_delimiters == 2){
                    if (token == "1"){
                        is_same_strand = true;
                    }
                    else if (token == "0"){
                        is_same_strand = false;
                    }
                    else{
                        throw runtime_error("ERROR: uninterpretable strand indicator on line ");
                    }

                    auto id_a = adjacency_map.insert_node(name_a);
                    auto id_b = adjacency_map.insert_node(name_b);

                    adjacency_map.edges.resize(adjacency_map.id_vs_name.left.size());

                    ShastaLabel label(true, false, false, false);
                    size_t index = adjacency_map.insert_edge(id_a, id_b, !is_same_strand, label);

                    if (alignment_files.size() <= index){
                        alignment_files.resize(index+1, "");
                    }

                    alignment_files[index] = file_iter.path();

                    // Stop iterating after the strand token
                    break;
                }

                token.resize(0);
                n_delimiters++;
            }
            else{
                token += c;
            }
        }
    }

    if (not (adjacency_map.labels.size() == alignment_files.size())){
        throw runtime_error("ERROR: adjacency size does not match alignment files: "
                            + to_string(adjacency_map.labels.size()) + " VS " + to_string(alignment_files.size()));
    }
}


void load_excluded_read_names_as_set(path excluded_reads_path, unordered_set<string>& excluded_reads) {
    ifstream file(excluded_reads_path);
    if (not file.good()){
        throw runtime_error("ERROR: excluded reads file could not be read: " + excluded_reads_path.string());
    }

    string line;

    while (getline(file, line)){
        excluded_reads.emplace(line.substr(0, line.size()));
    }
}


void evaluate_overlaps(
        path ref_overlap_path,
        path alignment_directory,
        path output_directory,
        uint32_t min_quality,
        path excluded_reads_path
){

    if (exists(output_directory)){
        throw runtime_error("ERROR: output directory already exists");
    }
    else{
        create_directories(output_directory);
    }

    AdjacencyMap adjacency;
    vector <path> alignment_files;

    load_alignment_files_as_adjacency_map(alignment_directory, adjacency, alignment_files);

    unordered_set<string> excluded_reads;

    if (not excluded_reads_path.empty()){
        // TODO: add chimera/palindrome labeling
        load_excluded_read_names_as_set(excluded_reads_path, excluded_reads);
    }

    for (const auto& name: excluded_reads){
        adjacency.erase_node(name);
    }

    add_paf_edges_to_adjacency_map(ref_overlap_path, min_quality, adjacency);

    path edge_csv_file_path = output_directory / "labeled_candidates.csv";
    write_edges_to_csv(edge_csv_file_path, adjacency);

    create_directories(output_directory / "true");
    create_directories(output_directory / "false");

    for (size_t i=0; i<adjacency.labels.size(); i++){
        // Ignore ref-only entries
        if (adjacency.labels[i].in_candidates){
            if (adjacency.labels[i].in_ref){
                cerr << '\n' << alignment_files[i] << '\n';
                path destination = output_directory / "true" / alignment_files[i].filename();
                rename(alignment_files[i], destination);
            }
            else{
                path destination = output_directory / "false" / alignment_files[i].filename();
                rename(alignment_files[i], destination);
            }
        }
    }
}


int main(int argc, char* argv[]){
    path ref_overlap_path;
    path alignment_directory;
    path output_directory;
    path excluded_reads_path;
    uint32_t min_quality;

    options_description options("Arguments:");

    options.add_options()
            ("a",
             value<path>(&ref_overlap_path)
                     ->required(),
             "REFERENCE GRAPH: path of file from which overlap can be inferred:\n"
             "\tPAF file containing alignments to some reference\n "
             "\tCSV file containing a list of reads pairs, with an indication for whether the overlap is cross-strand\n")

            ("alignment_directory",
             value<path>(&alignment_directory)
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
            alignment_directory,
            output_directory,
            min_quality,
            excluded_reads_path);

    return 0;


}

