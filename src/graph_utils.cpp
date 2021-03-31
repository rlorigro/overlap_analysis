#include "graph_utils.hpp"

namespace overlap_analysis {

void assign_graph_node_labels(
        Graph& graph,
        GraphAttributes& graph_attributes,
        vector<node>& nodes,
        uint32_string_bimap& id_vs_name,
        bool double_stranded,
        path output_path
) {
    ofstream file(output_path);

    graph_attributes.addAttributes(GraphAttributes::nodeLabel);

    for (size_t id = 0; id < nodes.size(); id++) {
        auto single_stranded_id = id;

        // If specified, undo the even/odd encoding scheme to reduce multiple node ids to a single mapping id->name
        if (double_stranded) {
            single_stranded_id = (id - (id % 2)) / 2;
        }

        if (nodes[id] == nullptr) {
            continue;
        }

        auto label = id_vs_name.left.at(single_stranded_id);

        // Keep track of whether this node was forward or positive and write it in the name
        if (double_stranded) {
            label += ((id % 2 > 0) ? "-" : "+");
        }

        if (output_path.empty()) {
            graph_attributes.label(nodes[id]) = label;
        }
            // Names might be too big to write in the graph, so optionally use numeric ids and write a separate csv table
        else {
            graph_attributes.label(nodes[id]) = to_string(id);

            file << id << ',' << label << '\n';
        }
    }
}


void get_all_read_names(
        Graph& graph,
        vector<node>& nodes,
        uint32_string_bimap& id_vs_name,
        bool double_stranded,
        bool strip_directional_tag,
        set<string>& read_names) {

    for (size_t id = 0; id < nodes.size(); id++) {
        auto single_stranded_id = id;

        // If specified, undo the even/odd encoding scheme to reduce multiple node ids to a single mapping id->name
        if (double_stranded) {
            single_stranded_id = (id - (id % 2)) / 2;
        }

        if (nodes[id] == nullptr) {
            continue;
        }

        auto label = id_vs_name.left.at(single_stranded_id);

        if (strip_directional_tag) {
            if (label.back() == '+' or label.back() == '-') {
                label = label.substr(0, label.size() - 1);
            }
        }

        read_names.emplace(label);
    }
}

}