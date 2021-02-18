#include "OverlapMap.hpp"
#include "PafElement.hpp"
#include "graph_utils.hpp"

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
#include <unordered_map>
#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <queue>

using std::experimental::filesystem::create_directories;
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
        ogdf::Graph& graph,
        vector<node>& nodes){

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
                            // Don't duplicate edges
                            if (graph.searchEdge(nodes[id], nodes[other_id]) == nullptr){
//                                cerr << "Creating edge: " << id << "->" << other_id << '\n';

                                uint32_t flipped_id;
                                uint32_t flipped_other_id;

                                if (id % 2 == 0){
                                    flipped_id = id + 1;
                                }
                                else{
                                    flipped_id = id - 1;
                                }

                                if (other_id % 2 == 0){
                                    flipped_other_id = other_id + 1;
                                }
                                else{
                                    flipped_other_id = other_id - 1;
                                }

                                graph.newEdge(nodes[id], nodes[other_id]);
                                graph.newEdge(nodes[flipped_id], nodes[flipped_other_id]);
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
/// -------------       ----                        ----
/// 1                   "Query sequence name"       string (assumed to be numeric for this project, using ONT reads)
/// 6                   "Target sequence name"      string
/// 8                   "Target start..."           int
/// 9                   "Target end..."             int
/// 12                  "Mapping quality"           int
///
void load_paf_as_graph(
        path paf_path,
        uint32_string_bimap& id_vs_name,
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
    string read_name;
    uint32_t start = 0;
    uint32_t stop = 0;
    uint32_t quality = 0;
    bool is_reverse = false;

    uint64_t n_delimiters = 0;
    uint64_t n_lines = 0;
    char c;

    while (paf_file.get(c)) {
        if (c == '\t') {
            if (n_delimiters == 0) {
                read_name = token;
            }
            else if (n_delimiters == 4) {
                is_reverse = (token == "-");
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

//                cerr << region_name << " " << start << " " << stop << '\n';

                if (quality >= min_quality) {
                    auto result = id_vs_name.right.find(read_name);
                    uint32_t id;
                    uint32_t forward_id;
                    uint32_t reverse_id;

                    // The bimap for id <-> name is not necessarily initialized for this read.
                    // If it does exist, then just fetch the ID
                    if (result != id_vs_name.right.end()) {
                        id = result->second;
                    }
                    // If it doesn't exist, then make a new ID by incrementing by 1 (aka get the size)
                    else{
                        id = id_vs_name.size();
                        id_vs_name.insert(bimap_pair(id, read_name));
                    }

                    forward_id = 2*id;
                    reverse_id = 2*id + 1;

                    if (not is_reverse) {
                        overlap_map.insert(region_name, start, stop, forward_id);
                    }
                    else {
                        overlap_map.insert(region_name, start, stop, reverse_id);
                    }

                    // In the case where a new read has just been encountered, add it to the end of the nodes list,
                    // which should maintain size = bimap.size()
                    if (reverse_id >= nodes.size()){
                        // Forward
                        nodes.emplace_back(graph.newNode());

                        // Reverse
                        nodes.emplace_back(graph.newNode());
                    }
                }
            }

            token.resize(0);
            n_delimiters++;
        }
        else if (c == '\n'){
            if (n_delimiters < 11){
                throw runtime_error("ERROR: file provided does not contain sufficient tab delimiters to be PAF");
            }

            token.resize(0);
            n_delimiters = 0;
            n_lines++;
        }
        else {
            token += c;
        }
    }

    create_graph_edges_from_overlap_map(overlap_map, graph, nodes);
}


void load_excluded_read_names_as_set(path excluded_reads_path, set<string>& excluded_reads){
    ifstream file(excluded_reads_path);
    string line;

    while(getline(file,line)){
        excluded_reads.emplace(line.substr(0,line.size()-1));
    }
}


void exclude_reads_from_graph(
        Graph& graph,
        vector<node>& nodes,
        path excluded_reads_path,
        uint32_string_bimap id_vs_name
){
    set<string> excluded_reads;
    load_excluded_read_names_as_set(excluded_reads_path, excluded_reads);

    for (size_t id=0; id<nodes.size(); id++){
        if (nodes[id] != nullptr){
            auto result = id_vs_name.left.find(id);

            if (result == id_vs_name.left.end()){
                continue;
            }

            if (excluded_reads.count(result->second) > 0) {
                graph.delNode(nodes[id]);
            }
        }
    }
}


void create_subgraph(
        const Graph& graph,
        Graph& subgraph,
        uint32_string_bimap& id_vs_name,
        uint32_string_bimap& subgraph_id_vs_name,
        vector<node>& nodes,
        vector<node>& subgraph_nodes,
        string start_name,
        uint32_t radius){

    cerr << "Extracting subgraph of radius " << radius << " around node " << start_name << '\n';

    uint32_t start_id;
    auto result = id_vs_name.right.find(start_name);

    if (result == id_vs_name.right.end()){
        throw runtime_error("ERROR: subgraph start node " + start_name + " not found in graph");
    }
    else{
        start_id = result->second;
    }

    uint32_t forward_start_id = 2*start_id;
    uint32_t reverse_start_id = 2*start_id + 1;

    queue <BfsQueueElement> q;
    unordered_map <size_t, size_t> distance;

    // Make new nodes for the subgraph
    auto forward_subgraph_start_node = subgraph.newNode();
    auto reverse_subgraph_start_node = subgraph.newNode();

    // Update the subgraph nodes list and Keep track of the name that the subgraph node should map to
    subgraph_id_vs_name.insert(bimap_pair(subgraph_nodes.size(), start_name + '+'));
    subgraph_nodes.emplace_back(forward_subgraph_start_node);

    subgraph_id_vs_name.insert(bimap_pair(subgraph_nodes.size(), start_name + '-'));
    subgraph_nodes.emplace_back(reverse_subgraph_start_node);

    // Initialize the queue and make a new node for the forward and reverse nodes
    q.emplace(nodes[forward_start_id], forward_subgraph_start_node);
    q.emplace(nodes[reverse_start_id], reverse_subgraph_start_node);

    // Update the distance map, which serves the double purpose of tracking which nodes are visited.
    // This map tracks distance in the original graph, using the original IDs.
    distance[forward_start_id] = 0;
    distance[reverse_start_id] = 0;

    while (not q.empty()){
        auto item = q.front();

        q.pop();

        for (auto adjacent_element: item.original_node->adjEntries){
            node other = adjacent_element->twinNode();
            auto d = distance.at(item.original_node->index());

            // If this node's distance would be <= radius, do the update
            if (d < radius) {
                uint32_t single_stranded_id = (other->index() - (other->index() % 2)) / 2;

                // Find the name that corresponded to this node in the original graph
                string name = id_vs_name.left.at(single_stranded_id);

                // Create a directionally labeled name, because these nodes will be randomly ordered
                // and their ID cant indicate strandedness
                string subgraph_name = name + ((other->index() % 2 > 0)? "-" : "+");

                // If it hasn't been visited before, create a new node and a new edge
                if (distance.count(other->index()) == 0) {

                    auto new_node = subgraph.newNode();
                    subgraph.newEdge(item.subgraph_node, new_node);

                    q.emplace(other, new_node);
                    distance[other->index()] = d + 1;

                    subgraph_id_vs_name.insert(bimap_pair(subgraph_nodes.size(), subgraph_name));
                    subgraph_nodes.emplace_back(new_node);
                }
                // If it has been visited, just create a new edge
                else{
                    uint32_t id = subgraph_id_vs_name.right.at(subgraph_name);

                    node a = item.subgraph_node;
                    node b = subgraph_nodes[id];

                    if (subgraph.searchEdge(a, b) == nullptr){
                        subgraph.newEdge(a, b);
                    }
                }
            }
        }
    }
}


void plot_graph(
        path paf_path,
        uint32_t min_quality,
        path excluded_reads_path,
        string subgraph_node_name,
        uint32_t subgraph_radius) {

    bool double_stranded_labeling = true;

    path output_directory = paf_path.parent_path();
    create_directories(output_directory);

    uint32_string_bimap id_vs_name;

    RegionalOverlapMap overlap_map;
    Graph overlap_graph;
    vector<node> nodes;

    load_paf_as_graph(paf_path, id_vs_name, overlap_map, overlap_graph, nodes, min_quality);

    if (not excluded_reads_path.empty()){
        exclude_reads_from_graph(overlap_graph, nodes, excluded_reads_path, id_vs_name);
    }

    if (not subgraph_node_name.empty()){
        Graph subgraph;
        uint32_string_bimap subgraph_id_vs_name;
        vector<node> subgraph_nodes;

        create_subgraph(
            overlap_graph,
            subgraph,
            id_vs_name,
            subgraph_id_vs_name,
            nodes,
            subgraph_nodes,
            subgraph_node_name,
            subgraph_radius);

        overlap_graph = subgraph;
        nodes = subgraph_nodes;
        id_vs_name = subgraph_id_vs_name;

        double_stranded_labeling = false;
    }

    GraphAttributes graph_attributes;
    assign_default_graph_rendering_attributes(overlap_graph, graph_attributes);
    assign_graph_node_labels(overlap_graph, graph_attributes, nodes, id_vs_name, double_stranded_labeling);

    write_graph_to_svg(overlap_graph, graph_attributes, paf_path.replace_extension("double_stranded.svg"));
}


pair<string, uint32_t> parse_subgraph_argument(string subgraph_argument){
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
    path excluded_reads_path;
    uint32_t min_quality;
    string subgraph_argument;
    string subgraph_node_name;
    uint32_t subgraph_radius;

    options_description options("Arguments:");

    options.add_options()
            ("paf_path",
             value<path>(&paf_path)
             ->required(),
             "File path of PAF file containing alignments to some reference")

            ("mapq",
            value<uint32_t>(&min_quality)
            ->required()
            ->default_value(50),
            "Minimum allowed mapping quality to load into the graph")

            ("exclude",
             value<path>(&excluded_reads_path)
             ->default_value(""),
             "File path of PAF file containing alignments to some reference")

            ("subgraph",
             value<string>(&subgraph_argument)
             ->default_value(""),
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

    if (not subgraph_argument.empty()) {
        tie(subgraph_node_name, subgraph_radius) = parse_subgraph_argument(subgraph_argument);
    }
    else{
        subgraph_node_name = "";
        subgraph_radius = 0;
    }

    plot_graph(paf_path, min_quality, excluded_reads_path, subgraph_node_name, subgraph_radius);

    return 0;
}
