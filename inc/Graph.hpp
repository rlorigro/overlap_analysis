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
#include <functional>
#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <queue>
#include <array>
#include <map>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::absolute;
using std::experimental::filesystem::exists;
using std::experimental::filesystem::path;
using std::runtime_error;
using std::unordered_map;
using std::function;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::string;
using std::vector;
using std::queue;
using std::pair;
using std::cerr;
using std::cout;
using std::array;
using std::map;

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


class ShastaLabel {
public:
    bool passes_readgraph2_criteria;
    bool in_read_graph;

    ShastaLabel()=default;
    ShastaLabel(bool passes_readgraph2_criteria, bool in_read_graph);
};


template <class T> class EdgeLabels {
public:
    vector <array <map <size_t, T>, 2> > data;
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
    virtual bool has_edge(const string& a, const string& b) const = 0;
    virtual bool has_node(const string& name) const = 0;

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

    bool has_node(const string& name) const;

    bool has_edge(const string& a, const string& b) const;

    // Not strictly needed except for methods that also accept double stranded graphs
    uint32_t get_single_stranded_id(uint32_t id) const;
};


class DoubleStrandedGraph: public UndirectedGraph {
public:
    EdgeLabels <ShastaLabel> edge_labels;

    /// Methods ///
    DoubleStrandedGraph()=default;

    uint32_t add_node(const string& read_name);
    void add_edge(uint32_t a, uint32_t b, bool allow_duplicates=true);
    bool remove_node(const string& name);

    uint32_t get_forward_id(uint32_t id) const;
    uint32_t get_reverse_id(uint32_t id) const;
    uint32_t get_single_stranded_id(uint32_t id) const;
    uint32_t get_double_stranded_id(uint32_t id, bool is_reverse) const;

    void get_all_read_names(set<string>& read_names);
    void assign_graph_node_labels(path output_path="");

    bool has_node(const string& name) const;

    bool has_edge(const string& a, const string& b) const;
    bool has_edge(const string& a, bool a_reversal, const string& b, bool b_reversal) const;

    bool is_reverse(uint32_t id) const;

    void node_union(const UndirectedGraph& other_graph);

    void create_subgraph(const string& start_name, uint32_t radius, Graph& subgraph);

    template <class T> void insert_label(uint32_t id0, uint32_t id1, bool is_cross_strand, T& label);
    template <class T> bool find_label(uint32_t id0, uint32_t id1, bool is_cross_strand, T& label);
};


class EdgeDescriptor {
public:
    const string& name0;
    const string& name1;
    const uint32_t id0;
    const uint32_t id1;
    const bool is_cross_strand;
    const bool in_ref;
    const bool in_non_ref;

    EdgeDescriptor(
        const string& name0,
        const string& name1,
        const uint32_t id0,
        const uint32_t id1,
        const bool is_cross_strand,
        const bool in_ref,
        const bool in_non_ref);
};


class EdgeDiff{
public:
    /// Attributes ///
    set <edge> a_only_edges;
    set <edge> b_only_edges;
    set <edge> a_both_edges;
    set <edge> b_both_edges;

    /// Methods ///
    EdgeDiff()=default;
    void agnostic_diff(const DoubleStrandedGraph& a, const DoubleStrandedGraph& b, path output_directory);
    void for_each_edge_comparison(
            const DoubleStrandedGraph& ref_graph,
            const DoubleStrandedGraph& graph,
            const function<void(EdgeDescriptor& e)>& f);
};


ostream& operator<<(ostream& o, EdgeDiff& g);


void create_graph_edges_from_overlap_map(
        RegionalOverlapMap& overlap_map,
        DoubleStrandedGraph& graph);


void load_paf_as_graph(
        path paf_path,
        RegionalOverlapMap& overlap_map,
        DoubleStrandedGraph& graph,
        uint32_t min_quality);


void load_adjacency_csv_as_graph(path adjacency_path, DoubleStrandedGraph& graph);


template <class T> void DoubleStrandedGraph::insert_label(uint32_t id0, uint32_t id1, bool is_cross_strand, T& label){
    if (id0 >= edge_labels.data.size()){
        edge_labels.data.resize(id0+1);
    }
    edge_labels.data[id0][is_cross_strand].emplace(id1,label);
}


template <class T> bool DoubleStrandedGraph::find_label(uint32_t id0, uint32_t id1, bool is_cross_strand, T& label){
    if (id0 >= edge_labels.data.size() or id1 >= edge_labels.data.size()){
        cerr << "WARNING: attempting to locate edge label with node ID > or == to number of nodes" << '\n';
        return false;
    }

    bool success = false;

    // Edges are bidirectional, but only one label is stored, so both directions must be searched
    auto iter0 = edge_labels.data[id0][is_cross_strand].find(id1);
    if (iter0 != edge_labels.data[id0][is_cross_strand].end()){
        label = iter0->second;
        success = true;
    }
    else{
        auto iter1 = edge_labels.data[id1][is_cross_strand].find(id0);
        if (iter1 != edge_labels.data[id1][is_cross_strand].end()){
            label = iter1->second;
            success = true;
        }
    }

    return success;
}


}

#endif //OVERLAP_ANALYSIS_GRAPH_H
