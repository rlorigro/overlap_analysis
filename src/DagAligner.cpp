#include "DiagonalTree.hpp"
#include "DagAligner.hpp"
#include "SvgPlot.hpp"

#include <numeric>

using std::max;


Node::Node(size_t x, size_t y, size_t length):
    prev(),
    next(),
    x(x),
    y(y),
    length(length),
    score(0)
{}


ostream& operator<<(ostream& o, Node& n){
    o << '(' << n.x << ',' << n.y << ')' << ' ' << n.length << ' ' << n.score << '\n';

    o << "Prev:" << '\n' << '\t';
    for (auto i: n.prev){
        o << i << ' ';
    }
    o << '\n';

    o << "Next:" << '\n' << '\t';
    for (auto i: n.next){
        o << i << ' ';
    }
    o << '\n';

    return o;
}


// TODO: Remove this code block?
Dag::Dag(const vector<xy_match>& matches, size_t x_size, size_t y_size, size_t max_gap):
        nodes(),
        alignment_path(),
        x_size(x_size),
        y_size(y_size),
        max_gap(max_gap)
{
    rtree <xy_match, quadratic<30> > tree(matches);

    for (auto& match: matches){
        size_t length = match.second;

        size_t x_start = match.first.get<0>();
        size_t x_stop = x_start + length;

        size_t y_start = match.first.get<1>();
        size_t y_stop = y_start + length;

        std::vector<xy_match> results;
        box q(xy_point{x_stop,y_stop},xy_point{x_stop+max_gap,y_stop+max_gap});

        tree.query(intersects(q), std::back_inserter(results));

        for (const auto& r: results){
            // TODO: add nodes to graph
            cerr << x_start << ' ' << y_start << " (" << r.first.get<0>() << ',' << r.first.get<1>() << ")=" << r.second << '\n';
        }
        cerr << '\n';
    }
}


Dag::Dag(const vector <pair <coord_t,size_t> >& matches, size_t x_size, size_t y_size, size_t max_gap):
        nodes(),
        alignment_path(),
        x_size(x_size),
        y_size(y_size),
        max_gap(max_gap)
{
    //TODO: add and test `first_only` flag

    DiagonalTree tree(x_size, y_size);

    // Node index 0 is reserved for the 'source' placeholder
    nodes.emplace_back(0,0,0);

    // Node index 1 is reserved for the 'sink' placeholder
    nodes.emplace_back(x_size,y_size,0);

    // Construct spatial index for matches in the matrix
    for (const auto& m: matches){
        auto x = m.first.first;
        auto y = m.first.second;

        nodes.emplace_back(x, y, m.second);

        // Each entry in the tree points to an index in the node list
        tree.insert(x, y, nodes.size() - 1);
    }

    // Construct graph based on nearest down-diagonal neighbors in the spatial index
    for (size_t i=0; i<nodes.size(); i++){
        if (i <= SINK_INDEX){
            continue;
        }

        auto& node = nodes[i];

        size_t x_stop = node.x + node.length;
        size_t y_stop = node.y + node.length;

        vector <pair <coord_t, size_t> > results;
        tree.find(x_stop, y_stop, max_gap, results);

        // Given the nearest down-diagonal matches within a certain gap size, add edges in the graph to those matches
        for (const auto& r: results){
            auto& adjacent_node = nodes[r.second];
            node.next.emplace_back(r.second);
            adjacent_node.prev.emplace_back(i);
            n_edges++;

            cerr << node.x << ' ' << node.y << " (" << r.first.first << ',' << r.first.second << ")=" << r.second << '\n';
        }
        cerr << '\n';

        // If it's close to the top/left boundary of the matrix then connect it to the source node
        if (node.x <= max_gap or node.y <= max_gap){
            nodes[SOURCE_INDEX].next.emplace_back(i);
            node.prev.emplace_back(SOURCE_INDEX);
            n_edges++;
        }

        // If it's close to the bottom/right boundary of the matrix then connect it to the sink node
        if ((x_stop >= x_size - max_gap) or (y_stop >= y_size - max_gap)){
            nodes[SINK_INDEX].prev.emplace_back(i);
            node.next.emplace_back(SINK_INDEX);
            n_edges++;
        }
    }
}


void Dag::compute_alignment(){
    vector <size_t> max_prev(nodes.size());

    // Forward pass
    for_each_in_kahns_iteration([&](const size_t i){
        int64_t max_score;

        if (i == SOURCE_INDEX){
            max_score = 0;
        }
        else {
            max_score = std::numeric_limits<int64_t>::min();

            for (const auto j: nodes[i].prev){
                if (nodes[j].score > max_score){
                    max_score = nodes[j].score;
                    max_prev[i] = j;
                }
            }
        }

        cerr << nodes[i].x << ' ' << nodes[i].y << ' ' << max_score << '\n';

    });

    for (size_t i=0; i<max_prev.size(); i++){
        cerr << i << ' ' << max_prev[i] << '\n';
    }

    // Traceback
    auto i = SINK_INDEX;
    while (i != SOURCE_INDEX){
        alignment_path.emplace_front(i);
        i = max_prev[i];
    }
}


void Dag::write_to_svg(path output_path){
    SvgPlot plot(output_path, 800, 800, 0, x_size, 0, y_size, true);
    string point_type = "circle";

    // Plot all the matches and their edges
    for (size_t i=0; i<nodes.size(); i++){
        auto& node = nodes[i];

        size_t x_stop = node.x + node.length;
        size_t y_stop = node.y + node.length;

        plot.add_line(node.x, node.y, x_stop, y_stop, double(max(x_size,y_size))/100, "black");
        plot.add_point(node.x, node.y, point_type, double(max(x_size,y_size))/80, "black");

        for (size_t j: node.prev){
            auto& prev_node = nodes[j];

            auto prev_x_stop = prev_node.x + prev_node.length;
            auto prev_y_stop = prev_node.y + prev_node.length;

            plot.add_line(prev_x_stop, prev_y_stop, node.x, node.y, double(max(x_size,y_size))/400, "orange");
        }
    }

    // Plot the alignment path that was computed if it exists
    for (size_t path_index=0; path_index<alignment_path.size(); path_index++){
        auto i = alignment_path[path_index];

        auto& node = nodes[i];

        size_t x_stop = node.x + node.length;
        size_t y_stop = node.y + node.length;

        plot.add_line(node.x, node.y, x_stop, y_stop, double(max(x_size,y_size))/600, "red");
        plot.add_point(node.x, node.y, point_type, double(max(x_size,y_size))/160, "red");

        if (path_index > 0){
            auto j = alignment_path[path_index-1];

            auto& prev_node = nodes[j];

            auto prev_x_stop = prev_node.x + prev_node.length;
            auto prev_y_stop = prev_node.y + prev_node.length;

            plot.add_line(prev_x_stop, prev_y_stop, node.x, node.y, double(max(x_size,y_size))/600, "red");
            plot.add_point(node.x, node.y, point_type, double(max(x_size,y_size))/160, "red");
        }
    }
}


ostream& operator<<(ostream& o, Dag& dag){

    size_t i = 0;
    for (auto& n: dag.nodes){
        o << "NODE " << i++ << '\n';
        o << n << '\n';
    }
    o << '\n';

    return o;
}
