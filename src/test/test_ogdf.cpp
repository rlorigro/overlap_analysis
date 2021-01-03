#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/fileformats/GraphIO.h>

using ogdf::Graph;
using ogdf::GraphAttributes;
using ogdf::node;
using ogdf::edge;
using ogdf::DPolyline;
using ogdf::DPoint;


/// Taken from https://ogdf.github.io/doc/ogdf/ex-basic.html
int main(){
    Graph G;
    GraphAttributes GA(G,GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics);

    const int LEN = 11;
    for(int i = 1; i < LEN; ++i) {
        node left = G.newNode();
        GA.x(left) = -5*(i+1);
        GA.y(left) = -20*i;
        GA.width(left) = 10*(i+1);
        GA.height(left) = 15;

        node bottom = G.newNode();
        GA.x(bottom) = 20*(LEN-i);
        GA.y(bottom) = 5*(LEN+1-i);
        GA.width(bottom) = 15;
        GA.height(bottom) = 10*(LEN+1-i);

        edge e = G.newEdge(left,bottom);
        DPolyline &p = GA.bends(e);
        p.pushBack(DPoint(10,-20*i));
        p.pushBack(DPoint(20*(LEN-i),-10));
    }

    ogdf::GraphIO::write(GA, "output-manual.gml", ogdf::GraphIO::writeGML);
    ogdf::GraphIO::write(GA, "output-manual.svg", ogdf::GraphIO::drawSVG);

    return 0;
}