#ifndef OVERLAP_ANALYSIS_GRAPH_UTILS_HPP
#define OVERLAP_ANALYSIS_GRAPH_UTILS_HPP

#include "OverlapMap.hpp"
#include "PafElement.hpp"

#include "boost/program_options.hpp"
#include "boost/bimap.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using boost::bimap;

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

#include <experimental/filesystem>
#include <iostream>
#include <utility>
#include <string>
#include <vector>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::path;
using std::runtime_error;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::string;
using std::vector;
using std::pair;
using std::cerr;
using std::cout;

typedef bimap<uint32_t,string> uint32_string_bimap;
typedef uint32_string_bimap::value_type bimap_pair;


class BfsQueueElement{
public:
    node original_node;
    node subgraph_node;

    BfsQueueElement()=default;
    BfsQueueElement(node original_node, node subgraph_node);
};


void assign_default_graph_rendering_attributes(Graph& graph, GraphAttributes& graph_attributes);


void assign_graph_node_labels(
        Graph& graph,
        GraphAttributes& graph_attributes,
        vector<node>& nodes,
        uint32_string_bimap& id_vs_name,
        bool double_stranded=false,
        path output_path={}
);

void write_graph_to_svg(Graph& graph, GraphAttributes& graph_attributes, path output_path);

void get_all_read_names(
        Graph& graph,
        vector<node>& nodes,
        uint32_string_bimap& id_vs_name,
        bool double_stranded,
        bool strip_directional_tag,
        set<string>& read_names);

#endif //OVERLAP_ANALYSIS_GRAPH_UTILS_HPP
