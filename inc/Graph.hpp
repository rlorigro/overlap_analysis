#ifndef OVERLAP_ANALYSIS_GRAPH_H
#define OVERLAP_ANALYSIS_GRAPH_H

#include "OverlapMap.hpp"
#include "PafElement.hpp"
#include "graph_utils.hpp"

#include "boost/program_options.hpp"
#include "boost/bimap.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::bool_switch;
using boost::program_options::value;
using boost::bimap;

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/graphics.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/fileformats/GraphIO.h>

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
using std::experimental::filesystem::absolute;
using std::experimental::filesystem::exists;
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


namespace overlap_analysis{


class BfsQueueElement{
public:
    node original_node;
    node subgraph_node;

    BfsQueueElement()=default;
    BfsQueueElement(node original_node, node subgraph_node);
};


class UndirectedGraph {
public:
    /// Attributes ///
    ogdf::Graph graph;
    vector<node> nodes;
    uint32_string_bimap id_vs_name;

    GraphAttributes graph_attributes;

    /// Methods ///
    virtual void assign_graph_node_labels(path output_path)=0;
    virtual bool has_edge(const string& a, const string& b) const=0;

    void assign_default_graph_rendering_attributes();
    void write_graph_to_svg(path output_path);

    virtual uint32_t get_single_stranded_id(uint32_t id) const=0;
};


class Graph: public UndirectedGraph{
public:
    /// Attributes ///
    bool pseudodirectional;

    /// Methods ///
    Graph()=default;
    void get_all_read_names(set<string>& read_names);
    void assign_graph_node_labels(path output_path="");

    bool has_edge(const string& a, const string& b) const;

    // Not strictly needed except for methods that also accept double stranded graphs
    uint32_t get_single_stranded_id(uint32_t id) const;
};


class DoubleStrandedGraph: public UndirectedGraph {
public:

    /// Methods ///
    DoubleStrandedGraph()=default;

    uint32_t add_node(const string& read_name);
    void add_edge(uint32_t a, uint32_t b);
    bool remove_node(const string& name);

    uint32_t get_forward_id(uint32_t id) const;
    uint32_t get_reverse_id(uint32_t id) const;
    uint32_t get_single_stranded_id(uint32_t id) const;
    uint32_t get_double_stranded_id(uint32_t id, bool is_reverse) const;

    void get_all_read_names(set<string>& read_names);
    void assign_graph_node_labels(path output_path="");

    bool has_edge(const string& a, const string& b) const;
    bool has_edge(const string& a, bool a_reversal, const string& b, bool b_reversal) const;

    bool is_reverse(uint32_t id) const;

    void create_subgraph(const string& start_name, uint32_t radius, Graph& subgraph);
};


class GraphDiff{
public:
    /// Attributes ///
    const UndirectedGraph& graph_a;
    const UndirectedGraph& graph_b;
    set <edge> a_only_edges;
    set <edge> b_only_edges;
    set <edge> a_both_edges;
    set <edge> b_both_edges;

    /// Methods ///
    GraphDiff(const DoubleStrandedGraph& a, const DoubleStrandedGraph& b);
};


ostream& operator<<(ostream& o, GraphDiff& g);


void create_graph_edges_from_overlap_map(
        RegionalOverlapMap& overlap_map,
        DoubleStrandedGraph& graph);


void load_paf_as_graph(
        path paf_path,
        RegionalOverlapMap& overlap_map,
        DoubleStrandedGraph& graph,
        uint32_t min_quality);


void load_adjacency_csv_as_graph(path adjacency_path, DoubleStrandedGraph& graph);


}

#endif //OVERLAP_ANALYSIS_GRAPH_H
