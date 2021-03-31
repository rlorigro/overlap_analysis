#include "Graph.hpp"


namespace overlap_analysis{


BfsQueueElement::BfsQueueElement(node original_node, node subgraph_node):
        original_node(original_node),
        subgraph_node(subgraph_node)
{};


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


uint32_t DoubleStrandedGraph::get_forward_id(uint32_t id){
    return 2*id;
}


uint32_t DoubleStrandedGraph::get_reverse_id(uint32_t id){
    return 2*id + 1;
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
        id_vs_name.insert(bimap_pair(id, read_name));
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


void DoubleStrandedGraph::add_edge(uint32_t a, uint32_t b){
    // Don't duplicate edges
    if (graph.searchEdge(nodes[a], nodes[b]) == nullptr){

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

        graph.newEdge(nodes[a], nodes[b]);
        graph.newEdge(nodes[flipped_id], nodes[flipped_other_id]);
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
    subgraph.id_vs_name.insert(bimap_pair(subgraph.nodes.size(), start_name + '+'));
    subgraph.nodes.emplace_back(forward_subgraph_start_node);

    subgraph.id_vs_name.insert(bimap_pair(subgraph.nodes.size(), start_name + '-'));
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

                    subgraph.id_vs_name.insert(bimap_pair(subgraph.nodes.size(), subgraph_name));
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
        auto single_stranded_id = id;

        // If specified, undo the even/odd encoding scheme to reduce multiple node ids to a single mapping id->name
        single_stranded_id = (id - (id % 2)) / 2;

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


}

