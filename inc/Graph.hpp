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
typedef bimap<size_t,string> sizet_string_bimap;


namespace overlap_analysis{


class Accuracy{
public:
    size_t n_false_positives = 0;
    size_t n_true_positives = 0;
    size_t n_false_negatives = 0;

    double compute_precision();
    double compute_sensitivity();
};


class BfsQueueElement{
public:
    node original_node;
    node subgraph_node;

    BfsQueueElement()=default;
    BfsQueueElement(node original_node, node subgraph_node);
};


class ShastaLabel {
public:
    bool in_candidates;
    bool passes_readgraph2_criteria;
    bool in_read_graph;
    bool in_ref;

    ShastaLabel();
    ShastaLabel(bool in_candidates, bool passes_readgraph2_criteria, bool in_read_graph, bool in_ref);

    ShastaLabel& operator|=(ShastaLabel& b){
        in_candidates = in_candidates or b.in_candidates;
        passes_readgraph2_criteria = passes_readgraph2_criteria or b.passes_readgraph2_criteria;
        in_read_graph = in_read_graph or b.in_read_graph;
        in_ref = in_ref or b.in_ref;

        return *this;
    }
};


class AdjacencyMap {
public:
    vector <array <map <size_t, size_t>, 2> > edges;
    sizet_string_bimap id_vs_name;
    vector <ShastaLabel> labels;

    size_t insert_node(const string& read_name);
    size_t insert_edge(uint32_t id0, uint32_t id1, bool is_cross_strand, ShastaLabel& label);
    bool find(uint32_t id0, uint32_t id1, bool is_cross_strand, ShastaLabel& label);
    bool find(string& name0, string& name1, bool is_cross_strand, ShastaLabel& label);
    pair<bool,size_t> find(uint32_t id0, uint32_t id1, bool is_cross_strand);
    pair<bool,size_t> find(string& name0, string& name1, bool is_cross_strand);
    void erase_edge(uint32_t id0, uint32_t id1, bool is_cross_strand);
    void erase_node(uint32_t id);
    void erase_node(const string& name);
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
};


void create_graph_edges_from_overlap_map(
        RegionalOverlapMap& overlap_map,
        DoubleStrandedGraph& graph);


void for_each_paf_element(
        path paf_path,
        uint32_t min_quality,
        const function<void(
                string& region_name,
                string& read_name,
                uint32_t start,
                uint32_t stop,
                uint32_t quality,
                bool is_reverse)>& f);


void for_each_paf_element(
        path paf_path,
        uint32_t min_quality,
        const function<void(PafElement& paf_element)>& f
);


void for_each_overlap_in_overlap_map(
        RegionalOverlapMap& overlap_map,
        const function<void(size_t id, size_t other_id)>& f
);


void load_paf_as_graph(
        path paf_path,
        DoubleStrandedGraph& graph,
        uint32_t min_quality);


void load_paf_as_adjacency_map(path paf_path, AdjacencyMap& adjacency_map, uint32_t min_quality);


void for_each_edge_in_shasta_adjacency_csv(
        path adjacency_path,
        const function<void(
                string& name_a,
                string& name_b,
                bool is_cross_strand,
                ShastaLabel& label)>& f);


void load_adjacency_csv_as_graph(path adjacency_path, DoubleStrandedGraph& graph);


void load_adjacency_csv_as_adjacency_map(path adjacency_path, AdjacencyMap& adjacency_map);


void add_paf_edges_to_adjacency_map(path paf_path, uint32_t min_quality, AdjacencyMap& adjacency_map);


void for_each_edge_in_adjacency(
        AdjacencyMap& a,
        const function<void(
                const string& name0,
                const string& name1,
                bool is_cross_strand,
                const ShastaLabel& label
        )>& f);


void write_edges_to_csv(path file_path, AdjacencyMap& a);

}

#endif //OVERLAP_ANALYSIS_GRAPH_H
