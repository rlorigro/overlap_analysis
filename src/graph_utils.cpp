#include "graph_utils.hpp"

BfsQueueElement::BfsQueueElement(node original_node, node subgraph_node):
    original_node(original_node),
    subgraph_node(subgraph_node)
{};


void assign_default_graph_rendering_attributes(Graph& graph, GraphAttributes& graph_attributes){
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


void assign_graph_node_labels(
        Graph& graph,
        GraphAttributes& graph_attributes,
        vector<node>& nodes,
        uint32_string_bimap& id_vs_name,
        bool double_stranded
){

    graph_attributes.addAttributes(GraphAttributes::nodeLabel);

    for (size_t id=0; id<nodes.size(); id++){
        auto single_stranded_id = id;

        if (double_stranded){
            single_stranded_id = (id - (id % 2)) / 2;
        }

        if (nodes[id] == nullptr){
            continue;
        }

        auto label = id_vs_name.left.at(single_stranded_id);

        if (double_stranded){
            label += ((id % 2 > 0)? "-" : "+");
        }

        graph_attributes.label(nodes[id]) = label;
    }
}


void write_graph_to_svg(Graph& graph, GraphAttributes& graph_attributes, path output_path){

    FMMMLayout layout_engine;
    layout_engine.useHighLevelOptions(true);
    layout_engine.unitEdgeLength(40.0);
    layout_engine.newInitialPlacement(true);
    layout_engine.qualityVersusSpeed(FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);

//    cerr << layout_engine.springStrength() << '\n';

    layout_engine.call(graph_attributes);

    graph_attributes.directed() = false;
    ogdf::GraphIO::write(graph_attributes, output_path, ogdf::GraphIO::drawSVG);

}
