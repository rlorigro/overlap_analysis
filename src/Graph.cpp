#include "OverlapMap.hpp"
#include "Graph.hpp"
#include "ogdf/basic/simple_graph_alg.h"

using ogdf::makeParallelFree;

namespace overlap_analysis{


ShastaLabel::ShastaLabel(bool passes_readgraph2_criteria, bool in_read_graph, bool in_ref):
    passes_readgraph2_criteria(passes_readgraph2_criteria),
    in_read_graph(in_read_graph),
    in_ref(in_ref)
{}


ShastaLabel::ShastaLabel():
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
        id = result->second;
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
        id = result->second;
    }
    // If it doesn't exist, then make a new ID by incrementing by 1 (aka get the size)
    else{
        id = id_vs_name.size();
        id_vs_name.insert(uint32_string_bimap::value_type(id, read_name));
    }

    reverse_id = 2*id + 1;

    // In the case where a new read has just been encountered, add it to the end of the nodes list,
    // which should maintain size = bimap.size()
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
        auto a_id = get_forward_id(a_iter->second);
        auto b_id = get_forward_id(b_iter->second);

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
        auto a_id = get_double_stranded_id(a_iter->second, a_reversal);
        auto b_id = get_double_stranded_id(b_iter->second, b_reversal);

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
    // Don't duplicate edges

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
        start_id = result->second;
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
void create_graph_edges_from_overlap_map(
        RegionalOverlapMap& overlap_map,
        DoubleStrandedGraph& graph) {

    set<uint32_t> empty_set = {};

    for (auto& item: overlap_map.intervals){
        const auto& overlaps = item.second;
        auto& prev_read_set = empty_set;

        for (auto i = begin(overlaps), e = end(overlaps); i!=e; ++i){
            auto& read_set = i->second;

            for (const auto& id: read_set){
                // If this read id is not in the previous set, it indicates that more edges need to be built
                if (prev_read_set.count(id) == 0){
                    for (const auto& other_id: read_set){
                        if (other_id != id){
                            graph.add_edge(id, other_id);
                        }
                    }
                }
            }

            prev_read_set = read_set;
        }
    }

    // Remove any duplicated edges
    makeParallelFree(graph.graph);
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

    cerr << "\tParsing PAF as interval map..." << '\n';
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

    cerr << "\tConstructing graph from interval map..." << '\n';
    create_graph_edges_from_overlap_map(overlap_map, graph);
}


void load_adjacency_csv_as_graph(path adjacency_path, DoubleStrandedGraph& graph){
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

                ShastaLabel l(passes_readgraph2_criteria, in_readgraph);
                graph.insert_label(id_a, id_b, is_cross_strand, l);
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


void load_adjacency_csv_as_adjacency_map(path adjacency_path, EdgeLabels<ShastaLabel>& adjacency_map){
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

                ShastaLabel l(passes_readgraph2_criteria, in_readgraph);

                auto id_a = adjacency_map.insert_node(name_a);
                auto id_b = adjacency_map.insert_node(name_b);
                adjacency_map.insert_edge(id_a, id_b, is_cross_strand, l);
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


EdgeDescriptor::EdgeDescriptor(
        const string& name0,
        const string& name1,
        const uint32_t id0,
        const uint32_t id1,
        const bool is_cross_strand,
        const bool in_ref,
        const bool in_non_ref):
    name0(name0),
    name1(name1),
    id0(id0),
    id1(id1),
    is_cross_strand(is_cross_strand),
    in_ref(in_ref),
    in_non_ref(in_non_ref)
{}


/// Dumb brute force search to compare edges in 2 graphs
/// A better method might be a joint BFS
void EdgeDiff::agnostic_diff(const DoubleStrandedGraph& a, const DoubleStrandedGraph& b, path output_directory){
    ofstream file_a(output_directory / "a_only_edges.csv");
    ofstream file_b(output_directory / "b_only_edges.csv");

    for (auto edge: a.graph.edges){
        auto nodes = edge->nodes();

        bool reversal0 = a.is_reverse(nodes[0]->index());
        bool reversal1 = a.is_reverse(nodes[1]->index());

        auto id0 = a.get_single_stranded_id(nodes[0]->index());
        auto id1 = a.get_single_stranded_id(nodes[1]->index());

        const auto& name0 = a.id_vs_name.left.at(id0);
        const auto& name1 = a.id_vs_name.left.at(id1);

        if (b.has_edge(name0, reversal0, name1, reversal1)){
            a_both_edges.insert(edge);
        }
        else{
            file_a << name0 << ',' << name1 << ',' << (reversal0 == reversal1) <<'\n';
            a_only_edges.insert(edge);
        }
    }

    for (auto edge: b.graph.edges){
        auto nodes = edge->nodes();

        bool reversal0 = b.is_reverse(nodes[0]->index());
        bool reversal1 = b.is_reverse(nodes[1]->index());

        auto id0 = b.get_single_stranded_id(nodes[0]->index());
        auto id1 = b.get_single_stranded_id(nodes[1]->index());

        const auto& name0 = b.id_vs_name.left.at(id0);
        const auto& name1 = b.id_vs_name.left.at(id1);

        if (a.has_edge(name0, reversal0, name1, reversal1)){
            b_both_edges.insert(edge);
        }
        else{
            file_b << name0 << ',' << name1 << ',' << (reversal0 == reversal1) <<'\n';
            b_only_edges.insert(edge);
        }
    }
}


/// Dumb brute force search to compare edges in 2 graphs
/// A better method might be a joint BFS
void EdgeDiff::for_each_edge_comparison(
        const DoubleStrandedGraph& ref_graph,
        const DoubleStrandedGraph& graph,
        const function<void(EdgeDescriptor& e)>& f
){
    for (auto edge: ref_graph.graph.edges){
        auto nodes = edge->nodes();

        bool reversal0 = ref_graph.is_reverse(nodes[0]->index());
        bool reversal1 = ref_graph.is_reverse(nodes[1]->index());

        auto id0 = ref_graph.get_single_stranded_id(nodes[0]->index());
        auto id1 = ref_graph.get_single_stranded_id(nodes[1]->index());

        const auto& name0 = ref_graph.id_vs_name.left.at(id0);
        const auto& name1 = ref_graph.id_vs_name.left.at(id1);

        bool is_cross_strand = (reversal0 != reversal1);
        bool in_ref = true;
        bool in_non_ref = graph.has_edge(name0, reversal0, name1, reversal1);

        EdgeDescriptor e(
            name0,
            name1,
            id0,
            id1,
            is_cross_strand,
            in_ref,
            in_non_ref);

        f(e);
    }

    for (auto edge: graph.graph.edges){
        auto nodes = edge->nodes();

        bool reversal0 = graph.is_reverse(nodes[0]->index());
        bool reversal1 = graph.is_reverse(nodes[1]->index());

        auto id0 = graph.get_single_stranded_id(nodes[0]->index());
        auto id1 = graph.get_single_stranded_id(nodes[1]->index());

        const auto& name0 = graph.id_vs_name.left.at(id0);
        const auto& name1 = graph.id_vs_name.left.at(id1);

        bool is_cross_strand = (reversal0 != reversal1);
        bool in_ref = ref_graph.has_edge(name0, reversal0, name1, reversal1);
        bool in_non_ref = true;

        // All the mutual edges have been iterated already, don't use them again
        if (not in_ref){
            EdgeDescriptor e(
                name0,
                name1,
                id0,
                id1,
                is_cross_strand,
                in_ref,
                in_non_ref);

            f(e);
        }
    }
}


ostream& operator<<(ostream& o, EdgeDiff& g){
    o << "a_edges," << g.a_only_edges.size() << '\n';
    o << "b_edges," << g.b_only_edges.size() << '\n';
    o << "a_and_b_edges," << g.b_both_edges.size() << '\n';

    return o;
}


}

