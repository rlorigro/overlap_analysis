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
using std::experimental::filesystem::exists;
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


namespace overlap_analysis {

template<class T> void render_graph(T& graph, uint16_t label_type, path output_directory, string name="graph") {
    path svg_path = output_directory / (name + ".svg");

    if (not exists(output_directory)){
        throw runtime_error("ERROR: output directory does not exist: " + output_directory.string());
    }

    cerr << "Assigning graph visual attributes...\n";
    graph.assign_default_graph_rendering_attributes();

    if (label_type == 1) {
        // Type 1 is to directly print the labels in the SVG
        graph.assign_graph_node_labels();
    } else if (label_type == 2) {
        // Type 2 prints only integer IDs and then creates a lookup table CSV separately for ID:name
        path label_csv_path = output_directory / "node_names.csv";
        graph.assign_graph_node_labels(label_csv_path);
    }

    cerr << "Computing graph layout and saving SVG\n";
    graph.write_graph_to_svg(svg_path);

}

}

#endif //OVERLAP_ANALYSIS_GRAPH_UTILS_HPP
