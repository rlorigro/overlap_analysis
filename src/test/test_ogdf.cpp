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

#include <stdexcept>
#include <vector>
#include <string>

using std::runtime_error;
using std::vector;
using std::string;


int main(){
    Graph g;
    GraphAttributes graph_attributes(g,
        GraphAttributes::nodeGraphics |
        GraphAttributes::edgeGraphics |
        GraphAttributes::nodeLabel |
        GraphAttributes::edgeStyle |
        GraphAttributes::nodeStyle |
        GraphAttributes::nodeTemplate);

    vector<node> nodes;

    size_t n_nodes = 5;

    for (size_t i=0; i<n_nodes; i++) {
        auto n = g.newNode(i);
        nodes.emplace_back(n);
        graph_attributes.shape(n) = Shape::Ellipse;
    }

    for (size_t i1=0; i1<n_nodes; i1++) {
        for (size_t i2=i1; i2<n_nodes; i2++) {
            if (i1 == i2) {
                continue;
            }

            g.newEdge(nodes[i1], nodes[i2]);
        }
    }

    FMMMLayout layout_engine;

    layout_engine.useHighLevelOptions(true);
    layout_engine.unitEdgeLength(15.0);
    layout_engine.newInitialPlacement(true);
    layout_engine.qualityVersusSpeed(FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);
    layout_engine.call(graph_attributes);

    graph_attributes.directed() = false;
    ogdf::GraphIO::write(graph_attributes, "test_ogdf.svg", ogdf::GraphIO::drawSVG);



    return 0;
}