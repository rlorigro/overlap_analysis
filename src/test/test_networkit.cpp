#include <networkit/graph/Graph.hpp>
#include <networkit/viz/MaxentStress.hpp>

using NetworKit::Graph;
using NetworKit::node;
using NetworKit::edgeid;
using NetworKit::MaxentStress;

#include <vector>
#include <string>
#include <stdexcept>

using std::vector;
using std::string;
using std::runtime_error;


void run_command(string& argument_string){
    int exit_code = system(argument_string.c_str());

    if (exit_code != 0){
        throw runtime_error("ERROR: command failed to run: " + argument_string);
    }
}


int main(){
    Graph g;

    vector<node> nodes;
    size_t n_nodes = 5;

    for (size_t i=0; i<n_nodes; i++){
        node n = g.addNode();
        nodes.emplace_back(n);
    }

    for (size_t i1=0; i1<n_nodes; i1++) {
        for (size_t i2 = i1; i2<n_nodes; i2++) {
            if (i1 == i2){
                continue;
            }
            g.addEdge(nodes[i1], nodes[i2]);
        }
    }

    size_t n_dimensions = 2;
    size_t k = 5;
    MaxentStress layout_engine(g, n_dimensions, k, 0.5);

    layout_engine.run();

    layout_engine.writeGraphToGML("test_networkit.gml");
    string command = "gml2gv test_networkit.gml | dot -Tsvg -o test_networkit.svg";
    run_command(command);

    return 0;
}
