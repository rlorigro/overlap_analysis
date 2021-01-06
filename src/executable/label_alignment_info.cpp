#include "OverlapMap.hpp"
#include "PafElement.hpp"

#include "boost/program_options.hpp"
#include "boost/bimap.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using boost::bimap;

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/graphics.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/energybased/SpringEmbedderKK.h>
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
void assign_graph_edges_from_overlap_map(
        RegionalOverlapMap& overlap_map,
        ogdf::Graph& graph,
        vector<node>& nodes){

    set<uint32_t> empty_set = {};

    for (auto& item: overlap_map.intervals){
        const auto& overlaps = item.second;
        auto& prev_read_set = empty_set;

        // TODO: Iterate the sets and add any edges that don't exist yet in the graph, as they are encountered.
        for (auto i = begin(overlaps), e = end(overlaps); i !=e; ++i){
            const auto& interval = i->first;
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
                            // Don't duplicate edges
                            if (graph.searchEdge(nodes[id], nodes[other_id]) == nullptr){
//                                cerr << "Creating edge: " << id << "->" << other_id << '\n';
                                graph.newEdge(nodes[id], nodes[other_id]);
                            }
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
        vector<node>& nodes,
        uint32_t min_quality){

    ifstream paf_file(paf_path);

    if (not paf_file.good()){
        throw runtime_error("ERROR: could not open input file: " + paf_path.string());
    }

    string token;
    string region_name;
    uint32_t read_name;     // Assume numeric read names
    uint32_t start;
    uint32_t stop;
    uint32_t quality;

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
            }
            else if (n_delimiters == 11) {
                quality = stoi(token);

                cerr << region_name << " " << start << " " << stop << '\n';

                if (quality >= min_quality) {
                    uint32_t id = id_vs_name.right.at(read_name);
                    overlap_map.insert(region_name, start, stop, id);

                    if (nodes[id] == nullptr) {
                        nodes[id] = graph.newNode();
                    }
                }
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

    assign_graph_edges_from_overlap_map(overlap_map, graph, nodes);
}


void write_graph_to_svg(Graph& graph, vector<node>& nodes, uint32_bimap& id_vs_name, path output_path, bool label=false){
    GraphAttributes graph_attributes(
            graph,
            GraphAttributes::nodeGraphics |
            GraphAttributes::edgeGraphics |
            GraphAttributes::nodeLabel |
            GraphAttributes::edgeStyle |
            GraphAttributes::nodeStyle |
            GraphAttributes::nodeTemplate);

    uint32_t node_diameter = 8;
    ogdf::Color edge_color(20, 20, 20, 255);

    for (size_t id=0; id<nodes.size(); id++){
        if (nodes[id] == nullptr){
            continue;
        }

        graph_attributes.shape(nodes[id]) = ogdf::Shape::Ellipse;
        graph_attributes.width(nodes[id]) = node_diameter;
        graph_attributes.height(nodes[id]) = node_diameter;

        if (label) {
            graph_attributes.label(nodes[id]) = to_string(id_vs_name.left.at(id));
        }
    }

    for (auto edge: graph.edges){
        graph_attributes.strokeColor(edge) = edge_color;
    }


    FMMMLayout layout_engine;
    layout_engine.useHighLevelOptions(true);
    layout_engine.unitEdgeLength(40.0);
    layout_engine.newInitialPlacement(true);
    layout_engine.qualityVersusSpeed(FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);

    cerr << layout_engine.springStrength() << '\n';

    layout_engine.call(graph_attributes);

    graph_attributes.directed() = false;
    ogdf::GraphIO::write(graph_attributes, output_path, ogdf::GraphIO::drawSVG);

}


void load_paf_as_table(
        path paf_path,
        uint32_bimap& id_vs_name,
        vector <vector <PafElement> >& paf_table){

    ifstream paf_file(paf_path);

    if (not paf_file.good()){
        throw runtime_error("ERROR: could not open input file: " + paf_path.string());
    }

    string token;
    string region_name;
    uint32_t read_name;     // Assume numeric read names
    uint32_t start;
    uint32_t stop;
    uint32_t quality;

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
            }
            else if (n_delimiters == 11) {
                quality = stoi(token);

                uint32_t id = id_vs_name.right.at(read_name);
                paf_table[id].emplace_back(region_name, start, stop, quality);
                cerr << region_name << " " << start << " " << stop << '\n';
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


uint32_t compute_overlap(
        uint32_t id1,
        uint32_t id2,
        Graph& overlap_graph,
        vector<node>& nodes,
        vector <vector <PafElement> >& paf_table){

    uint32_t overlap_size = 0;

    if ((nodes[id1] != nullptr) and (nodes[id2] != nullptr)) {
        auto edge = overlap_graph.searchEdge(nodes[id1], nodes[id2]);

        if (edge != nullptr) {
            RegionalOverlapMap overlap;

            // Add all intervals from read 1
            for (auto& paf_element: paf_table[id1]) {
                overlap.insert(paf_element.target_name, paf_element.start, paf_element.stop, id1);
            }
            // Add all intervals from read 2
            for (auto& paf_element: paf_table[id2]) {
                overlap.insert(paf_element.target_name, paf_element.start, paf_element.stop, id2);
            }

            for (auto& region: overlap.intervals) {
                for (auto& item: region.second) {
                    auto& interval = item.first;
                    auto& ids = item.second;

                    if (ids.size() == 2) {
                        overlap_size += interval.upper() - interval.lower();
                    }
                }
            }
        }
    }

    return overlap_size;
}


/// \param info_csv_path
/// \param output_path
///
/// Parse the AlignmentInfo CSV with columns:
/// readId0,readId1,minAlignedFraction,markerCount,maxDrift,maxSkip,trim
///
/// and append a column which contains the overlap length between reads
///
void write_labelled_alignment_info(
        path info_csv_path,
        path output_path,
        Graph& overlap_graph,
        vector<node>& nodes,
        vector <vector <PafElement> >& paf_table){

    ifstream paf_file(info_csv_path);

    if (not paf_file.good()){
        throw runtime_error("ERROR: could not open input file: " + info_csv_path.string());
    }

    ofstream out_file(output_path);

    if (not out_file.is_open()){
        throw runtime_error("ERROR: could not write to output file: " + output_path.string());
    }


    string line;
    string token;
    string region_name;
    uint32_t read_id1;     // Assume numeric read names
    uint32_t read_id2;     // Assume numeric read names

    uint64_t n_delimiters = 0;
    uint64_t n_lines = 0;
    char c;

    while (paf_file.get(c)) {
        if (c == ',') {
            if (n_lines == 0){
                continue;
            }

            if (n_delimiters == 0) {
                read_id1 = stoi(token);
            }
            if (n_delimiters == 1) {
                read_id2 = stoi(token);
            }

            token.resize(0);
            n_delimiters++;
        }
        else if (c == '\n'){
            // Write header or append the line
            if (n_lines == 0){
                out_file << line << ",overlap";
            }
            else {
                uint32_t overlap_size = compute_overlap(read_id1, read_id2, overlap_graph, nodes, paf_table);
                out_file << line << "," << overlap_size;
            }

            token.resize(0);
            line.resize(0);
            n_delimiters = 0;
            n_lines++;
        }
        else {
            // Update the token if not a delimiter or newline char (i.e. not whitespace)
            token += c;
        }

        // Dump every character into the line string
        line += c;
    }
}


void label_alignment_info(path info_csv_path, path read_csv_path, path paf_path, path output_path){
    path output_directory = output_path.parent_path();
    create_directories(output_directory);

    uint32_bimap id_vs_name;
    uint32_t max_id;
    max_id = load_csv_as_id_map(read_csv_path, id_vs_name);

    RegionalOverlapMap overlap_map;
    Graph overlap_graph;
    vector<node> nodes;
    nodes.resize(max_id+1, nullptr);

    uint32_t min_quality = 50;
    load_paf_as_graph(paf_path, id_vs_name, overlap_map, overlap_graph, nodes, min_quality);

    vector <vector <PafElement> > paf_table;
    paf_table.resize(max_id+1);
    load_paf_as_table(paf_path, id_vs_name, paf_table);

    write_labelled_alignment_info(info_csv_path, output_path, overlap_graph, nodes, paf_table);

    write_graph_to_svg(overlap_graph, nodes, id_vs_name, "test_label_alignment_info.svg");
}


int main(int argc, char* argv[]){
    path info_csv_path;
    path read_csv_path;
    path paf_path;
    path output_path;

    options_description options("Arguments:");

    options.add_options()
            ("info_csv_path",
             value<path>(&info_csv_path)->required(),
             "File path of CSV alignment info dump created by shasta readGraph.creationMethod 2 with debug mode on")

            ("read_csv_path",
             value<path>(&read_csv_path)->required(),
             "File path of CSV file ReadSummary.csv produced by shasta")

            ("paf_path",
             value<path>(&paf_path)->required(),
             "File path of PAF file containing alignments to some reference")

            ("output_path",
             value<path>(&output_path)->required(),
             "File path dictating where to save the labeled alignment info CSV. Directories will be created.")
            ;

    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cerr << options << "\n";
        return 0;
    }
    notify(vm);

    label_alignment_info(info_csv_path, read_csv_path, paf_path, output_path);

    return 0;
}
