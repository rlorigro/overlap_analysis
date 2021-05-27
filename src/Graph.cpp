#include "OverlapMap.hpp"
#include "Graph.hpp"
#include "ogdf/basic/simple_graph_alg.h"

using ogdf::makeParallelFree;
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
void for_each_overlap_in_overlap_map(
        RegionalOverlapMap& overlap_map,
        const function<void(size_t id, size_t other_id)>& f
        ) {

    for (auto& item: overlap_map.intervals){
        cerr << item.first << '\n';
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

                if (paf_element.map_quality >= min_quality){
                    f(paf_element);
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
        id = result->second;
    }
    // If it doesn't exist, then make a new ID by incrementing by 1 (aka get the size)
    else{
        id = id_vs_name.size();
        id_vs_name.insert(sizet_string_bimap::value_type(id, read_name));
    }

    return id;
}


void AdjacencyMap::insert_edge(uint32_t id0, uint32_t id1, bool is_cross_strand, ShastaLabel& label){
    auto iter0 = edges.at(id0)[is_cross_strand].find(id1);

    // If an entry is not found, it should not be found in both directions
    if (iter0 == edges[id0][is_cross_strand].end()){
        labels.emplace_back(label);
        edges.at(id0)[is_cross_strand].emplace(id1,labels.size()-1);
        edges.at(id1)[is_cross_strand].emplace(id0,labels.size()-1);
    }
    // If the label/edge exists already, increment it (find the union of label membership)
    else if (iter0 != edges.at(id0)[is_cross_strand].end()){
        labels.at(iter0->second) |= label;
    }
    else{
        throw runtime_error("ERROR: asymmetrical entry in undirected graph " + to_string(id0) + " " + to_string(id1) + (is_cross_strand ? "0" : "1"));
    }
}


void AdjacencyMap::erase_edge(uint32_t id0, uint32_t id1, bool is_cross_strand){
    auto iter0 = edges.at(id0)[is_cross_strand].find(id1);

    // If an entry is found, erase both edge directions
    if (iter0 != edges[id0][is_cross_strand].end()){
        edges.at(id0)[is_cross_strand].erase(iter0);
        edges.at(id1)[is_cross_strand].erase(id0);
    }
}


void AdjacencyMap::erase_node(uint32_t id){
    // Find all edges associated with this node and remove them
    for (auto is_cross_strand: {0,1}){
        for (auto& item: edges[id][is_cross_strand]){
            auto id1 = item.first;

            erase_edge(id1,id,is_cross_strand);
        }
    }

    // Remove the name and id from the bimap
    id_vs_name.left.erase(id);
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
    bool success = this->find(iter0->second, iter1->second, is_cross_strand, label);

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

    for (size_t id0=0; id0<a.edges.size(); id0++) {
        string name0 = a.id_vs_name.left.at(id0);

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

}

