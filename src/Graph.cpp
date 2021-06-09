#include "OverlapMap.hpp"
#include "Graph.hpp"
#include "SvgPlot.hpp"
#include "ogdf/basic/simple_graph_alg.h"

using ogdf::makeParallelFree;

#include <map>

using std::multimap;
using std::tie;

namespace overlap_analysis{


ShastaLabel::ShastaLabel(bool in_candidates, bool passes_readgraph2_criteria, bool in_read_graph, bool in_ref):
    in_candidates(in_candidates),
    passes_readgraph2_criteria(passes_readgraph2_criteria),
    in_read_graph(in_read_graph),
    in_ref(in_ref)
{}


ShastaLabel::ShastaLabel():
    in_candidates(false),
    passes_readgraph2_criteria(false),
    in_read_graph(false),
    in_ref(false)
{}


BfsQueueElement::BfsQueueElement(node original_node, node subgraph_node):
    original_node(original_node),
    subgraph_node(subgraph_node)
{}


bool DoubleStrandedGraph::remove_node(const string& name){
    size_t id;
    auto result = id_vs_name.right.find(name);

    if (result == id_vs_name.right.end()){
        // Exit early and indicate that the node did not exist
        return false;
    }
    else{
        id = result->get_left();
    }

    // ID map is single stranded, but the nodes vector has 2 nodes for every 1 read
    uint32_t forward_id = 2*id;
    uint32_t reverse_id = 2*id + 1;

    if (nodes[forward_id] != nullptr) {
        graph.delNode(nodes[forward_id]);
        nodes[forward_id] = nullptr;
    }
    if (nodes[reverse_id] != nullptr) {
        graph.delNode(nodes[reverse_id]);
        nodes[reverse_id] = nullptr;
    }

    return true;
}


uint32_t DoubleStrandedGraph::get_forward_id(uint32_t id) const{
    return 2*id;
}


uint32_t DoubleStrandedGraph::get_reverse_id(uint32_t id) const{
    return 2*id + 1;
}


uint32_t Graph::get_single_stranded_id(uint32_t id) const{
    return id;
}


uint32_t DoubleStrandedGraph::get_single_stranded_id(uint32_t id) const{
    return (id - (id % 2)) / 2;
}


bool DoubleStrandedGraph::is_reverse(uint32_t id) const{
    return ((id % 2) == 1);
}


uint32_t DoubleStrandedGraph::get_double_stranded_id(uint32_t id, bool is_reverse) const{
    return 2*id + is_reverse;
}


uint32_t DoubleStrandedGraph::add_node(const string& read_name){
    auto result = id_vs_name.right.find(read_name);
    uint32_t id;
    uint32_t reverse_id;

    // The bimap for id <-> name is not necessarily initialized for this read.
    // If it does exist, then just fetch the ID
    if (result != id_vs_name.right.end()) {
        id = result->get_left();
    }
    // If it doesn't exist, then make a new ID by incrementing by 1 (aka get the size)
    else{
        id = id_vs_name.size();
        id_vs_name.insert(uint32_string_bimap::value_type(id, read_name));
    }

    reverse_id = 2*id + 1;

    // In the case where a new read has just been encountered, add it to the end of the nodes list,
    // which should maintain size = 2*bimap.size()
    if (reverse_id >= nodes.size()){
        // Forward
        nodes.emplace_back(graph.newNode());

        // Reverse
        nodes.emplace_back(graph.newNode());
    }

    return id;
}


bool DoubleStrandedGraph::has_edge(const string& a, const string& b) const{
    bool result = true;

    auto a_iter = id_vs_name.right.find(a);
    auto b_iter = id_vs_name.right.find(b);

    if (a_iter == id_vs_name.right.end() or b_iter == id_vs_name.right.end()){
        result = false;
    }
    else {
        auto a_id = get_forward_id(a_iter->get_left());
        auto b_id = get_forward_id(b_iter->get_left());

        if (graph.searchEdge(nodes[a_id], nodes[b_id]) == nullptr) {
            result = false;
        }
    }

    return result;
}


bool DoubleStrandedGraph::has_node(const string& name) const{
    bool result = false;

    auto iter = id_vs_name.right.find(name);

    if (iter != id_vs_name.right.end()){
        result = true;
    }

    return result;
}


bool Graph::has_node(const string& name) const{
    bool result = false;

    auto iter = id_vs_name.right.find(name);

    if (iter != id_vs_name.right.end()){
        result = true;
    }

    return result;
}


bool DoubleStrandedGraph::has_edge(const string& a, bool a_reversal, const string& b, bool b_reversal) const{
    bool result = true;

    auto a_iter = id_vs_name.right.find(a);
    auto b_iter = id_vs_name.right.find(b);

    if (a_iter == id_vs_name.right.end() or b_iter == id_vs_name.right.end()){
        result = false;
    }
    else {
        auto a_id = get_double_stranded_id(a_iter->get_left(), a_reversal);
        auto b_id = get_double_stranded_id(b_iter->get_left(), b_reversal);

        if (graph.searchEdge(nodes[a_id], nodes[b_id]) == nullptr) {
            result = false;
        }
    }

    return result;
}


bool Graph::has_edge(const string& a, const string& b) const{
    bool result = true;

    auto a_iter = id_vs_name.right.find(a);
    auto b_iter = id_vs_name.right.find(b);

    if (a_iter == id_vs_name.right.end() or b_iter == id_vs_name.right.end()){
        result = false;
    }
    else {
        auto a_id = a_iter->second;
        auto b_id = b_iter->second;

        if (graph.searchEdge(nodes[a_id], nodes[b_id]) == nullptr) {
            result = false;
        }
    }

    return result;
}


void DoubleStrandedGraph::add_edge(uint32_t a, uint32_t b, bool allow_duplicates){
    uint32_t flipped_id;
    uint32_t flipped_other_id;

    if (a % 2 == 0){
        flipped_id = a + 1;
    }
    else{
        flipped_id = a - 1;
    }

    if (b % 2 == 0){
        flipped_other_id = b + 1;
    }
    else{
        flipped_other_id = b - 1;
    }

    if (allow_duplicates) {
        graph.newEdge(nodes[a], nodes[b]);
        graph.newEdge(nodes[flipped_id], nodes[flipped_other_id]);
    }
    else{
        if (graph.searchEdge(nodes[a], nodes[b]) == nullptr){
            graph.newEdge(nodes[a], nodes[b]);
        }

        if (graph.searchEdge(nodes[flipped_id], nodes[flipped_other_id]) == nullptr){
            graph.newEdge(nodes[flipped_id], nodes[flipped_other_id]);
        }
    }
}


/// This method iterates outward from a source node using BFS to generate a single stranded subgraph
/// However it is single stranded only in the sense that the mapping from id -> node is no longer
/// id % 2 + int(is_reverse). Node names become read_name with a (+/-) suffix
void DoubleStrandedGraph::create_subgraph(const string& start_name, uint32_t radius, Graph& subgraph){
    // Need to remember how names were created so that +/- suffix can be stripped
    subgraph.pseudodirectional = true;

    uint32_t start_id;
    auto result = id_vs_name.right.find(start_name);

    if (result == id_vs_name.right.end()){
        throw runtime_error("ERROR: subgraph start node " + start_name + " not found in graph");
    }
    else{
        start_id = result->get_left();
    }

    uint32_t forward_start_id = 2*start_id;
    uint32_t reverse_start_id = 2*start_id + 1;

    queue <BfsQueueElement> q;
    unordered_map <size_t, size_t> distance;

    // Make new nodes for the subgraph
    auto forward_subgraph_start_node = subgraph.graph.newNode();
    auto reverse_subgraph_start_node = subgraph.graph.newNode();

    // Update the subgraph nodes list and Keep track of the name that the subgraph node should map to
    subgraph.id_vs_name.insert(uint32_string_bimap::value_type(subgraph.nodes.size(), start_name + '+'));
    subgraph.nodes.emplace_back(forward_subgraph_start_node);

    subgraph.id_vs_name.insert(uint32_string_bimap::value_type(subgraph.nodes.size(), start_name + '-'));
    subgraph.nodes.emplace_back(reverse_subgraph_start_node);

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

                    auto new_node = subgraph.graph.newNode();
                    subgraph.graph.newEdge(item.subgraph_node, new_node);

                    q.emplace(other, new_node);
                    distance[other->index()] = d + 1;

                    subgraph.id_vs_name.insert(uint32_string_bimap::value_type(subgraph.nodes.size(), subgraph_name));
                    subgraph.nodes.emplace_back(new_node);
                }
                // If it has been visited, just create a new edge
                else{
                    uint32_t id = subgraph.id_vs_name.right.at(subgraph_name);

                    node a = item.subgraph_node;
                    node b = subgraph.nodes[id];

                    if (subgraph.graph.searchEdge(a, b) == nullptr){
                        subgraph.graph.newEdge(a, b);
                    }
                }
            }
        }
    }
}


void Graph::get_all_read_names(set<string>& read_names){

    for (size_t id=0; id<nodes.size(); id++){
        auto single_stranded_id = id;

        if (nodes[id] == nullptr){
            continue;
        }

        auto label = id_vs_name.left.at(single_stranded_id);

        if (pseudodirectional){
            if (label.back() == '+' or label.back() == '-'){
                label = label.substr(0,label.size()-1);
            }
        }

        read_names.emplace(label);
    }
}


void DoubleStrandedGraph::get_all_read_names(set<string>& read_names){

    for (size_t id=0; id<nodes.size(); id++){
        auto single_stranded_id = id;

        // Undo the even/odd encoding scheme to reduce multiple node ids to a single mapping id->name
        single_stranded_id = (id - (id % 2)) / 2;

        if (nodes[id] == nullptr){
            continue;
        }

        auto label = id_vs_name.left.at(single_stranded_id);

        read_names.emplace(label);
    }
}


void DoubleStrandedGraph::node_union(const UndirectedGraph& other_graph){
    for (size_t id=0; id<nodes.size(); id++){
        auto single_stranded_id = id;

        // Undo the even/odd encoding scheme to reduce multiple node ids to a single mapping id->name
        single_stranded_id = (id - (id % 2)) / 2;

        if (nodes[id] == nullptr){
            continue;
        }

        auto label = id_vs_name.left.at(single_stranded_id);

        if (not other_graph.has_node(label)){
            graph.delNode(nodes[id]);
            nodes[id] = nullptr;
        }
    }
}


void UndirectedGraph::assign_default_graph_rendering_attributes(){
    graph_attributes = GraphAttributes(
            graph,
            GraphAttributes::nodeGraphics |
            GraphAttributes::edgeGraphics |
            GraphAttributes::edgeStyle |
            GraphAttributes::nodeStyle |
            GraphAttributes::nodeTemplate);

    uint32_t node_diameter = 8;
    ogdf::Color edge_color(20, 20, 20, 255);

    for (auto node: graph.nodes){
        if (node == nullptr){
            continue;
        }

        graph_attributes.shape(node) = ogdf::Shape::Ellipse;
        graph_attributes.width(node) = node_diameter;
        graph_attributes.height(node) = node_diameter;
    }

    for (auto edge: graph.edges){
        graph_attributes.strokeColor(edge) = edge_color;
        graph_attributes.strokeWidth(edge) = 0.3;
    }
}


void UndirectedGraph::write_graph_to_svg(path output_path) {

    FMMMLayout layout_engine;
    layout_engine.useHighLevelOptions(true);
    layout_engine.unitEdgeLength(40.0);
    layout_engine.newInitialPlacement(true);
    layout_engine.qualityVersusSpeed(FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);

//    cerr << layout_engine.springStrength() << '\n';

    layout_engine.call(graph_attributes);

    graph_attributes.directed() = false;

    cerr << "Writing svg: " << output_path << '\n';

    ogdf::GraphIO::write(graph_attributes, output_path, ogdf::GraphIO::drawSVG);
}


void DoubleStrandedGraph::assign_graph_node_labels(path output_path){
    ofstream file;

    if (not output_path.empty()){
        file.open(output_path);

        if (not file.good()){
            throw runtime_error("ERROR: could not write to file: " + output_path.string());
        }
    }

    graph_attributes.addAttributes(GraphAttributes::nodeLabel);

    for (size_t id=0; id<nodes.size(); id++){
        // Undo the even/odd encoding scheme to reduce multiple node ids to a single mapping id->name
        auto single_stranded_id = (id - (id % 2)) / 2;

        if (nodes[id] == nullptr){
            continue;
        }

        auto label = id_vs_name.left.at(single_stranded_id);

        // Keep track of whether this node was forward or positive and write it in the name
        label += ((id % 2 > 0)? "-" : "+");

        if (output_path.empty()) {
            graph_attributes.label(nodes[id]) = label;
        }
        // Names might be too big to write in the graph, so optionally use numeric ids and write a separate csv table
        else{
            graph_attributes.label(nodes[id]) = to_string(id);

            file << id << ',' << label << '\n';
        }
    }
}


void Graph::assign_graph_node_labels(path output_path){
    ofstream file;

    if (not output_path.empty()){
        file.open(output_path);

        if (not file.good()){
            throw runtime_error("ERROR: could not write to file: " + output_path.string());
        }
    }

    graph_attributes.addAttributes(GraphAttributes::nodeLabel);

    for (size_t id=0; id<nodes.size(); id++){
        auto single_stranded_id = id;

        if (nodes[id] == nullptr){
            continue;
        }

        auto label = id_vs_name.left.at(single_stranded_id);

        if (output_path.empty()) {
            graph_attributes.label(nodes[id]) = label;
        }
        // Names might be too big to write in the graph, so optionally use numeric ids and write a separate csv table
        else{
            graph_attributes.label(nodes[id]) = to_string(id);

            file << id << ',' << label << '\n';
        }
    }
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
void for_each_overlap_in_overlap_map(
        RegionalOverlapMap& overlap_map,
        const function<void(size_t id, size_t other_id)>& f
        ) {

    for (auto& item: overlap_map.intervals){
        const auto& overlaps = item.second;
        set<uint32_t> prev_read_set = {};

        for (auto i = begin(overlaps), e = end(overlaps); i!=e; ++i){
            auto& read_set = i->second;

            for (const auto& id: read_set){
                // If this read id is not in the previous set, it indicates that more edges need to be built
                if (prev_read_set.count(id) == 0){
                    for (const auto& other_id: read_set){
                        if (other_id != id){
                            f(id, other_id);
                        }
                    }
                }
            }

            prev_read_set = read_set;
        }
    }
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
void for_each_paf_element(
        path paf_path,
        uint32_t min_quality,
        const function<void(PafElement& paf_element)>& f
) {

    ifstream paf_file(paf_path);

    if (not paf_file.good()) {
        throw runtime_error("ERROR: could not open input file: " + paf_path.string());
    }

    string token;
    PafElement paf_element;

    uint64_t n_delimiters = 0;
    uint64_t n_lines = 0;
    char c;

    while (paf_file.get(c)){
        if (c == '\t') {
            if (n_delimiters == 0){
                paf_element.query_name = token;
            }
            else if (n_delimiters == 4){
                paf_element.is_reverse = (token == "-");
            }
            else if (n_delimiters == 5){
                paf_element.target_name = token;
            }
            else if (n_delimiters == 7){
                paf_element.start = stoi(token);
            }
            else if (n_delimiters == 8){
                paf_element.stop = stoi(token);
            }
            else if (n_delimiters == 11){
                paf_element.map_quality = stoi(token);
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

            if (paf_element.map_quality >= min_quality){
                f(paf_element);
            }

            token.resize(0);
            n_delimiters = 0;
            n_lines++;
        }
        else {
            token += c;
        }
    }
}


void load_paf_as_graph(
        path paf_path,
        DoubleStrandedGraph& graph,
        uint32_t min_quality) {

    RegionalOverlapMap overlap_map;

    for_each_paf_element(paf_path, min_quality, [&](PafElement& paf_element){
        uint32_t id;
        uint32_t forward_id;
        uint32_t reverse_id;

        id = graph.add_node(paf_element.query_name);

        forward_id = 2 * id;
        reverse_id = 2 * id + 1;

        if (not paf_element.is_reverse) {
            overlap_map.insert(paf_element.target_name, paf_element.start, paf_element.stop, forward_id);
        } else {
            overlap_map.insert(paf_element.target_name, paf_element.start, paf_element.stop, reverse_id);
        }
    });

    cerr << "\tConstructing graph from interval map..." << '\n';
    for_each_overlap_in_overlap_map(overlap_map, [&](size_t id, size_t other_id){
        graph.add_edge(id, other_id);
    });

    // Remove any duplicated edges
    makeParallelFree(graph.graph);
}


void add_paf_edges_to_adjacency_map(path paf_path, uint32_t min_quality, AdjacencyMap& adjacency_map){
    RegionalOverlapMap overlap_map;

    cerr << "\tParsing PAF as interval map..." << '\n';
    for_each_paf_element(paf_path, min_quality, [&](PafElement& paf_element){
        uint32_t id;
        uint32_t forward_id;
        uint32_t reverse_id;

        auto iter = adjacency_map.id_vs_name.right.find(paf_element.query_name);

        // Skip any entries for which the node doesn't exist already
        if (iter == adjacency_map.id_vs_name.right.end()){
//            cerr << "SKIPPING: " << paf_element.query_name << " because not found in graph\n";
            return;
        }

        id = iter->get_left();

        forward_id = 2 * id;
        reverse_id = 2 * id + 1;

//        cerr << "Found ref entry for " << paf_element.query_name << ' ' << id << ' ' << forward_id << ' ' << reverse_id << '\n';

        if (not paf_element.is_reverse) {
            overlap_map.insert(paf_element.target_name, paf_element.start, paf_element.stop, forward_id);
        } else {
            overlap_map.insert(paf_element.target_name, paf_element.start, paf_element.stop, reverse_id);
        }
    });

    cerr << "\tUpdating graph edges with inferred overlaps..." << '\n';

    for_each_overlap_in_overlap_map(overlap_map, [&](size_t id, size_t other_id){
        // Infer cross-strandedness using the even/odd id encoding
        bool is_cross_strand = ((id % 2) != (other_id % 2));

        // Convert back to single-stranded id
        id = (id - (id % 2)) / 2;
        other_id = (other_id - (other_id % 2)) / 2;

        // Create a label that only indicates ref membership
        ShastaLabel label(false, false, false, true);

        // Paf intervals can overlap (and will not be reduced in the set of ids if it is cross-strand)
        if (id == other_id){
            return;
        }

        // Insert or update the edge with ref membership
        adjacency_map.insert_edge(id, other_id, is_cross_strand, label);
    });
}


void load_paf_as_adjacency_map(path paf_path, AdjacencyMap& adjacency_map, uint32_t min_quality){

    RegionalOverlapMap overlap_map;

    for_each_paf_element(paf_path, min_quality, [&](PafElement& paf_element){
        uint32_t id;
        uint32_t forward_id;
        uint32_t reverse_id;

        id = adjacency_map.insert_node(paf_element.query_name);

        forward_id = 2 * id;
        reverse_id = 2 * id + 1;

        if (not paf_element.is_reverse) {
            overlap_map.insert(paf_element.target_name, paf_element.start, paf_element.stop, forward_id);
        } else {
            overlap_map.insert(paf_element.target_name, paf_element.start, paf_element.stop, reverse_id);
        }
    });

    cerr << "\tConstructing graph from interval map..." << '\n';

    adjacency_map.edges.resize(adjacency_map.id_vs_name.left.size());

    for_each_overlap_in_overlap_map(overlap_map, [&](size_t id, size_t other_id){
        // Infer cross-strandedness using the even/odd id encoding
        bool is_cross_strand = ((id % 2) != (other_id % 2));

        // Convert back to single-stranded id
        id = (id - (id % 2)) / 2;
        other_id = (other_id - (other_id % 2)) / 2;

        // Create a label that only indicates ref membership
        ShastaLabel label(false, false, false, true);

        adjacency_map.insert_edge(id, other_id, is_cross_strand, label);
    });
}


void load_adjacency_csv_as_graph(path adjacency_path, DoubleStrandedGraph& graph){
    for_each_edge_in_shasta_adjacency_csv(adjacency_path, [&](string& name_a,
                                                              string& name_b,
                                                              bool is_cross_strand,
                                                              ShastaLabel& label){
        auto id_a = graph.add_node(name_a);
        auto id_b = graph.add_node(name_b);

        if (not is_cross_strand) {
            graph.add_edge(graph.get_forward_id(id_a), graph.get_forward_id(id_b));
            graph.add_edge(graph.get_reverse_id(id_a), graph.get_reverse_id(id_b));
        }
        else{
            graph.add_edge(graph.get_forward_id(id_a), graph.get_reverse_id(id_b));
            graph.add_edge(graph.get_reverse_id(id_a), graph.get_forward_id(id_b));
        }
    });
}


void load_adjacency_csv_as_adjacency_map(path adjacency_path, AdjacencyMap& adjacency_map){
    for_each_edge_in_shasta_adjacency_csv(adjacency_path, [&](string& name_a,
                                                              string& name_b,
                                                              bool is_cross_strand,
                                                              ShastaLabel& label){
        auto id_a = adjacency_map.insert_node(name_a);
        auto id_b = adjacency_map.insert_node(name_b);

        adjacency_map.edges.resize(adjacency_map.id_vs_name.left.size());

        adjacency_map.insert_edge(id_a, id_b, is_cross_strand, label);
    });
}


void for_each_edge_in_shasta_adjacency_csv(
        path adjacency_path,
        const function<void(
                string& name_a,
                string& name_b,
                bool is_cross_strand,
                ShastaLabel& label)>& f
){
    ifstream file(adjacency_path);

    if (not file.good()){
        throw runtime_error("ERROR: could not read file: " + adjacency_path.string());
    }

    uint32_t n_delimiters = 0;
    uint32_t n_lines = 0;
    string token;
    string name_a;
    string name_b;
    bool is_cross_strand;
    bool passes_readgraph2_criteria = false;
    bool in_readgraph = false;

    char c;

    while (file.get(c)){
        if (c == ',') {
            if (n_delimiters == 0){
                name_a = token;
            }
            else if (n_delimiters == 1){
                name_b = token;
            }
            else if (n_delimiters == 2){
                // Depending on the format there may be more data after this token (shasta uses 'isSameStrand'=Yes|No)
                is_cross_strand = (token == "No");
            }
            else if (n_delimiters == 3){
                // Depending on the format there may be more data after this token
                passes_readgraph2_criteria = (token == "Yes");
            }
            else if (n_delimiters == 4){
                // Depending on the format there may be more data after this token
                in_readgraph = (token == "Yes");
            }

            token.resize(0);
            n_delimiters++;
        }
        else if (c == '\n'){
            // Skip header line
            if (n_lines != 0){
                if (n_delimiters < 2){
                    throw runtime_error(
                            "ERROR: file provided does not contain sufficient delimiters to be adjacency csv at line: " +
                            to_string(n_lines));
                }
                if (n_delimiters == 2){
                    // shasta uses 'isSameStrand'=Yes|No
                    is_cross_strand = (token == "No");
                }

                ShastaLabel l(true, passes_readgraph2_criteria, in_readgraph, false);

                f(name_a, name_b, is_cross_strand, l);
            }

            token.resize(0);
            n_delimiters = 0;
            n_lines++;
        }
        else {
            token += c;
        }
    }

}


size_t AdjacencyMap::insert_node(const string& read_name){
    auto result = id_vs_name.right.find(read_name);
    size_t id;

    // The bimap for id <-> name is not necessarily initialized for this read.
    // If it already exists, then just fetch the ID
    if (result != id_vs_name.right.end()) {
        id = result->get_left();
    }
    // If it doesn't exist, then make a new ID by incrementing by 1 (aka get the size)
    else{
        id = id_vs_name.size();
        id_vs_name.insert(sizet_string_bimap::value_type(id, read_name));
    }

    return id;
}


size_t AdjacencyMap::insert_edge(uint32_t id0, uint32_t id1, bool is_cross_strand, ShastaLabel& label){
    if (id0 == id1){
        throw runtime_error("ERROR: self edges not allowed: " + to_string(id0) + "->" + to_string(id1));
    }

    // Find the edge if it exists. The iterator points to a size_t label_index that contains the label for
    // this edge in labels[index]
    auto iter0 = edges.at(id0).at(is_cross_strand).find(id1);
    size_t edge_label_index;

    // If an entry is not found, it should not be found in both directions
    if (iter0 == edges[id0][is_cross_strand].end()){
        labels.emplace_back(label);
        edges.at(id0)[is_cross_strand].emplace(id1,labels.size()-1);
        edges.at(id1)[is_cross_strand].emplace(id0,labels.size()-1);

        // The edge label index is now the last index of the label vector
        edge_label_index = labels.size() - 1;
    }
    // If the label/edge exists already, increment it (find the union of label membership)
    else if (iter0 != edges.at(id0)[is_cross_strand].end()){
        labels.at(iter0->second) |= label;
        edge_label_index = iter0->second;
    }
    else{
        throw runtime_error("ERROR: asymmetrical entry in undirected graph " + to_string(id0) + " " + to_string(id1) + (is_cross_strand ? "0" : "1"));
    }

    return edge_label_index;
}


void AdjacencyMap::erase_edge(uint32_t id0, uint32_t id1, bool is_cross_strand){
    auto iter0 = edges.at(id0)[is_cross_strand].find(id1);

    // If an entry is found, erase both edge directions
    if (iter0 != edges[id0][is_cross_strand].end()){
        labels[iter0->second] = ShastaLabel(false,false,false,false);

        edges.at(id0)[is_cross_strand].erase(iter0);
        edges.at(id1)[is_cross_strand].erase(id0);
    }
}


void AdjacencyMap::erase_node(uint32_t id){
    // Find all edges associated with this node and remove them
    for (auto is_cross_strand: {0,1}){

        // Find all the other nodes that point to this id and delete those references
        for (auto& item: edges[id][is_cross_strand]){
            auto id1 = item.first;
            auto label_index = item.second;

            edges[id1][is_cross_strand].erase(id);

            // Re-label the edge so it is not a member of any class
            labels[label_index] = ShastaLabel(false,false,false,false);
        }
    }

    // Delete the entire set of outgoing edges from this node
    edges[id][false] = {};
    edges[id][true] = {};

    // Remove the name and id from the bimap
    id_vs_name.left.erase(id);
}


void AdjacencyMap::erase_node(const string& name){
    auto iter = id_vs_name.right.find(name);

    if (iter != id_vs_name.right.end()){
        erase_node(iter->get_left());
    }
}


bool AdjacencyMap::find(uint32_t id0, uint32_t id1, bool is_cross_strand, ShastaLabel& label) {
    if (id0 >= edges.size() or id1 >= edges.size()){
        cerr << "WARNING: attempting to locate edge label with node ID > or == to number of nodes" << '\n';
        return false;
    }

    bool success = false;

    auto iter0 = edges[id0][is_cross_strand].find(id1);
    if (iter0 != edges[id0][is_cross_strand].end()){
        label = labels.at(iter0->second);
        success = true;
    }

    return success;
}


pair<bool,size_t> AdjacencyMap::find(uint32_t id0, uint32_t id1, bool is_cross_strand) {
    if (id0 >= edges.size() or id1 >= edges.size()){
        cerr << "WARNING: attempting to locate edge label with node ID > or == to number of nodes" << '\n';
        return {false,0};
    }

    pair<bool,size_t> result = {false,0};

    auto iter0 = edges[id0][is_cross_strand].find(id1);
    if (iter0 != edges[id0][is_cross_strand].end()){
        result.first = true;
        result.second = iter0->second;
    }

    return result;
}


bool AdjacencyMap::find(string& name0, string& name1, bool is_cross_strand, ShastaLabel& label) {
    // First check if name0 exists
    auto iter0 = id_vs_name.right.find(name0);

    if (iter0 == id_vs_name.right.end()){
        return false;
    }

    // Then check if name1 exists
    auto iter1 = id_vs_name.right.find(name1);

    if (iter1 == id_vs_name.right.end()){
        return false;
    }

    // See if the edge between name0,id0 and name1,id1 exists
    bool success = this->find(iter0->get_left(), iter1->get_left(), is_cross_strand, label);

    return success;
}


pair<bool,size_t> AdjacencyMap::find(string& name0, string& name1, bool is_cross_strand) {
    // First check if name0 exists
    auto iter0 = id_vs_name.right.find(name0);

    if (iter0 == id_vs_name.right.end()){
        return {false,0};
    }

    // Then check if name1 exists
    auto iter1 = id_vs_name.right.find(name1);

    if (iter1 == id_vs_name.right.end()){
        return {false,0};
    }

    // See if the edge between name0,id0 and name1,id1 exists
    auto result = this->find(iter0->second, iter1->second, is_cross_strand);

    return result;
}


pair<double,double> get_forward_node_coordinate(DoubleStrandedGraph& graph, uint32_t id){
    auto forward_id = graph.get_forward_id(id);

    auto forward_node = graph.nodes.at(forward_id);

    double x = 0;
    double y = 0;

    if (forward_node != nullptr){
        x = graph.graph_attributes.x(forward_node);
        y = graph.graph_attributes.y(forward_node);
    }

    return {x,y};
}


pair<double,double> get_reverse_node_coordinate(DoubleStrandedGraph& graph, uint32_t id){
    auto reverse_id = graph.get_reverse_id(id);
    auto reverse_node = graph.nodes[reverse_id];

    double x = 0;
    double y = 0;

    if (reverse_node != nullptr){
        x = graph.graph_attributes.x(reverse_node);
        y = graph.graph_attributes.y(reverse_node);
    }

    return {x,y};
}


/*
 * COLORS:
 *   Red        #FF2800
 *   Orange     #FF9E00
 *   Yellow     #FFDD00
 *   Lime       #B8F400
 *   Green      #00C442
 *   Blue       #0658C2
 *   Indigo     #450BBA
 *   Violet     #C50094
 */

void get_edge_color(const ShastaLabel& label, string& color){
    if (label.in_read_graph){
        if (label.in_ref){
            color = "#00C442";      // Green
        }
        else {
            color = "#FF2800";      // Red
        }
    }
    else if (label.passes_readgraph2_criteria){
        if (label.in_ref){
            color = "#0658C2";      // Blue
        }
        else {
            color = "#FF9E00";      // Orange
        }
    }
    else if (label.in_candidates){
        if (label.in_ref){
            color = "#450BBA";      // Indigo
        }
        else {
            color = "#FFDD00";      // Yellow
        }
    }
    else if (label.in_ref){
        color = "#C50094";          // Violet
    }
}


size_t get_edge_z_order(const ShastaLabel& label){
    size_t z_order;

    if (label.in_read_graph){
        if (label.in_ref){
            z_order = 3;
        }
    }
    else if (label.passes_readgraph2_criteria){
        if (label.in_ref){
            z_order = 2;
        }
    }
    else if (label.in_candidates){
        if (label.in_ref){
            z_order = 1;
        }
    }
    else if (label.in_ref){
        z_order = 4;
    }

    return z_order;
}


class EdgeObject{
public:
    double x0;
    double y0;
    double x1;
    double y1;
    string color;

    EdgeObject()=default;

    EdgeObject(
            double x0,
            double y0,
            double x1,
            double y1,
            string& color
    );

};


EdgeObject::EdgeObject(
        double x0,
        double y0,
        double x1,
        double y1,
        string& color
):
        x0(x0),
        y0(y0),
        x1(x1),
        y1(y1),
        color(color)
{}



void AdjacencyMap::plot(path output_path, bool read_graph_only){
    DoubleStrandedGraph graph;

    for_each_edge_in_adjacency(*this, [&](const string& name0,
                                          const string& name1,
                                          bool is_cross_strand,
                                          const ShastaLabel& label){

        if (read_graph_only and not label.in_read_graph){
            return;
        }

        auto id_a = graph.add_node(name0);
        auto id_b = graph.add_node(name1);

//        cerr << id_a << " " << name0 << " " << graph.get_forward_id(id_a) << " " << graph.get_reverse_id(id_a) << '\n';
//        cerr << id_b << " " << name1 << " " << graph.get_forward_id(id_b) << " " << graph.get_reverse_id(id_b) << '\n';

        if (not is_cross_strand) {
            graph.add_edge(graph.get_forward_id(id_a), graph.get_forward_id(id_b));
            graph.add_edge(graph.get_reverse_id(id_a), graph.get_reverse_id(id_b));
        }
        else{
            graph.add_edge(graph.get_forward_id(id_a), graph.get_reverse_id(id_b));
            graph.add_edge(graph.get_reverse_id(id_a), graph.get_forward_id(id_b));
        }
    });

    graph.graph_attributes = GraphAttributes(graph.graph);

    FMMMLayout layout_engine;
    layout_engine.useHighLevelOptions(true);
    layout_engine.unitEdgeLength(1.0);
    layout_engine.newInitialPlacement(true);
    layout_engine.qualityVersusSpeed(FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);

    layout_engine.call(graph.graph_attributes);

    graph.graph_attributes.directed() = false;

    double x_max = 0;
    double y_max = 0;

    for (auto& [id, name]: graph.id_vs_name){
        double forward_coord_x0;
        double forward_coord_y0;

        double reverse_coord_x0;
        double reverse_coord_y0;

        tie(forward_coord_x0,forward_coord_y0) = get_forward_node_coordinate(graph,id);
        tie(reverse_coord_x0,reverse_coord_y0) = get_reverse_node_coordinate(graph,id);

        if (forward_coord_x0 > x_max){
            x_max = forward_coord_x0;
        }
        if (forward_coord_y0 > y_max){
            y_max = forward_coord_y0;
        }
        if (reverse_coord_x0 > x_max){
            x_max = reverse_coord_x0;
        }
        if (reverse_coord_y0 > y_max){
            y_max = reverse_coord_y0;
        }
    }

    SvgPlot plot(output_path, 1200, 1200, 0, x_max, 0, y_max);

    string edge_color;
    string node_color = "black";
    string type = "circle";

    cerr << graph.nodes.size() << " " << graph.id_vs_name.size() << '\n';

    map <size_t,vector<EdgeObject> > edge_objects;

    for_each_edge_in_adjacency(*this, [&](const string& name0,
                                          const string& name1,
                                          bool is_cross_strand,
                                          const ShastaLabel& label){

        if (read_graph_only and not label.in_read_graph){
            return;
        }

        auto id0 = graph.id_vs_name.right.at(name0);
        auto id1 = graph.id_vs_name.right.at(name1);

//        cerr << id0 << "___" << name0 << " " << id1 << "___" << name1 << " " << is_cross_strand << '\n';

        double forward_coord_x0;
        double forward_coord_y0;

        double reverse_coord_x0;
        double reverse_coord_y0;

        double forward_coord_x1;
        double forward_coord_y1;

        double reverse_coord_x1;
        double reverse_coord_y1;

        tie(forward_coord_x0,forward_coord_y0) = get_forward_node_coordinate(graph,id0);
        tie(reverse_coord_x0,reverse_coord_y0) = get_reverse_node_coordinate(graph,id0);

        tie(forward_coord_x1,forward_coord_y1) = get_forward_node_coordinate(graph,id1);
        tie(reverse_coord_x1,reverse_coord_y1) = get_reverse_node_coordinate(graph,id1);

        get_edge_color(label, edge_color);
        auto z_order = get_edge_z_order(label);

//        cerr << z_order << ' '
//             << edge_color << ' '
//             << label.in_candidates << ' '
//             << label.passes_readgraph2_criteria << ' '
//             << label.in_read_graph << ' '
//             << label.in_ref << '\n';

        EdgeObject e0;
        EdgeObject e1;

        if (is_cross_strand){
            e0 = EdgeObject(forward_coord_x0, forward_coord_y0, reverse_coord_x1, reverse_coord_y1, edge_color);
            e1 = EdgeObject(reverse_coord_x0, reverse_coord_y0, forward_coord_x1, forward_coord_y1, edge_color);
        }
        else{
            e0 = EdgeObject(forward_coord_x0, forward_coord_y0, forward_coord_x1, forward_coord_y1, edge_color);
            e1 = EdgeObject(reverse_coord_x0, reverse_coord_y0, reverse_coord_x1, reverse_coord_y1, edge_color);
        }

        edge_objects[z_order].emplace_back(e0);
        edge_objects[z_order].emplace_back(e1);
    });

    // Add the edges to the SVG in the order that is determined by their ShastaLabel
    for (auto& [z,item]: edge_objects){
        for (auto& e: item){
            plot.add_line(e.x0, e.y0, e.x1, e.y1, 0.5, e.color);
        }
    }

    for (auto& [id, name]: graph.id_vs_name){
        double forward_coord_x0;
        double forward_coord_y0;

        double reverse_coord_x0;
        double reverse_coord_y0;

        tie(forward_coord_x0,forward_coord_y0) = get_forward_node_coordinate(graph,id);
        tie(reverse_coord_x0,reverse_coord_y0) = get_reverse_node_coordinate(graph,id);

        string title_forward = name + "+";
        string title_reverse = name + "-";

        plot.add_point(forward_coord_x0, forward_coord_y0, type, 1, node_color, title_forward);
        plot.add_point(reverse_coord_x0, reverse_coord_y0, type, 1, node_color, title_reverse);
    }
}


void for_each_edge_in_adjacency(
        AdjacencyMap& a,
        const function<void(
                const string& name0,
                const string& name1,
                bool is_cross_strand,
                const ShastaLabel& label
        )>& f){

    // Initialize vector to track whether bidirectional edges have been visited
    vector <bool> visited_a(a.labels.size(), false);

    for (auto& [id0,name0]: a.id_vs_name) {
        for (bool is_cross_strand: {0, 1}) {
            for (auto& iter: a.edges[id0][is_cross_strand]) {
                auto id1 = iter.first;
                auto label_index_a = iter.second;

                // If this edge was visited already in A, then skip it
                if (visited_a.at(label_index_a)){
                    continue;
                }

                // Mark the edge as visited in this graph
                visited_a.at(label_index_a) = true;

                string name1 = a.id_vs_name.left.at(id1);

                f(name0, name1, is_cross_strand, a.labels[label_index_a]);
            }
        }
    }
}


double Accuracy::compute_precision(){
    return (double(n_true_positives) / double(n_true_positives + n_false_positives));
}


double Accuracy::compute_sensitivity(){
    return (double(n_true_positives) / double(n_true_positives + n_false_negatives));
}


void write_edges_to_csv(path file_path, AdjacencyMap& a){
    ofstream file(file_path);
    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: file could not written: " + file_path.string());
    }

    file << "name0,name1,is_same_strand,in_candidates,passes_readgraph2_criteria,in_read_graph,in_ref," << '\n';

    Accuracy candidate_accuracy;
    Accuracy filter_criteria_accuracy;
    Accuracy read_graph_accuracy;

    for_each_edge_in_adjacency(a,
                               [&](const string& name0,
                                   const string& name1,
                                   bool is_cross_strand,
                                   const ShastaLabel& label){
       file << name0 << ',' << name1 << ','
            << (is_cross_strand ? "No" : "Yes") << ','
            << (label.in_candidates ? "Yes" : "No") << ','
            << (label.passes_readgraph2_criteria ? "Yes" : "No") << ','
            << (label.in_read_graph ? "Yes" : "No") << ','
            << (label.in_ref ? "Yes" : "No") << ',' << '\n';

       if (label.in_candidates){
           if (label.in_ref){
               candidate_accuracy.n_true_positives++;
           }
           else {
               candidate_accuracy.n_false_positives++;
           }
       }
       if (label.passes_readgraph2_criteria){
           if (label.in_ref){
               filter_criteria_accuracy.n_true_positives++;
           }
           else {
               filter_criteria_accuracy.n_false_positives++;
           }
       }
       if (label.in_read_graph){
           if (label.in_ref){
               read_graph_accuracy.n_true_positives++;
           }
           else {
               read_graph_accuracy.n_false_positives++;
           }
       }
       if (label.in_ref){
           if (not label.in_candidates){
               candidate_accuracy.n_false_negatives++;
           }
           if (not label.passes_readgraph2_criteria){
               filter_criteria_accuracy.n_false_negatives++;
           }
           if (not label.in_read_graph){
               read_graph_accuracy.n_false_negatives++;
           }
       }
    });

    cerr << "candidate_accuracy:\n"
         << "true_positives: \t" << candidate_accuracy.n_true_positives << '\n'
         << "false_positives:\t" << candidate_accuracy.n_false_positives << '\n'
         << "false_negatives:\t" << candidate_accuracy.n_false_negatives << '\n'
         << "precision:\t" << candidate_accuracy.compute_precision() << '\n'
         << "sensitivity:\t" << candidate_accuracy.compute_sensitivity() << '\n';

    cerr << '\n';
    cerr << "filter_criteria_accuracy:\n"
         << "true_positives: \t" << filter_criteria_accuracy.n_true_positives << '\n'
         << "false_positives:\t" << filter_criteria_accuracy.n_false_positives << '\n'
         << "false_negatives:\t" << filter_criteria_accuracy.n_false_negatives << '\n'
         << "precision:\t" << filter_criteria_accuracy.compute_precision() << '\n'
         << "sensitivity:\t" << filter_criteria_accuracy.compute_sensitivity() << '\n';

    cerr << '\n';
    cerr << "read_graph_accuracy:\n"
         << "true_positives: \t" << read_graph_accuracy.n_true_positives << '\n'
         << "false_positives:\t" << read_graph_accuracy.n_false_positives << '\n'
         << "false_negatives:\t" << read_graph_accuracy.n_false_negatives << '\n'
         << "precision:\t" << read_graph_accuracy.compute_precision() << '\n'
         << "sensitivity:\t" << read_graph_accuracy.compute_sensitivity() << '\n';
}



}

