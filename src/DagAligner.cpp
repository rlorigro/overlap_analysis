#include "DiagonalIntervalTree.hpp"
#include "DiagonalTree.hpp"
#include "DagAligner.hpp"
#include "SvgPlot.hpp"

#include <cmath>

using std::max;


Node::Node(size_t x, size_t y, size_t length, size_t score):
    prev(),
    next(),
    x(x),
    y(y),
    length(length),
    score(score)
{}


size_t Node::get_x_stop() {
    if (length == 0){
        return x;
    }
    else {
        return x + length - 1;
    }
}


size_t Node::get_y_stop() {
    if (length == 0){
        return y;
    }
    else {
        return y + length - 1;
    }
}


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
        n_edges(0),
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
        n_edges(0),
        max_gap(max_gap)
{
    // Nodes are initialized with the largest possible negative score
    int64_t start_score = gap_open_score + x_size*gap_score + y_size*gap_score;

    DiagonalTree tree(x_size, y_size);

    // Node index 0 is reserved for the 'source' placeholder
    nodes.emplace_back(0,0,0,0);

    // Node index 1 is reserved for the 'sink' placeholder
    nodes.emplace_back(x_size,y_size,0,start_score);

    // Construct spatial index for matches in the matrix
    for (const auto& m: matches){
        auto x = m.first.first;
        auto y = m.first.second;

        nodes.emplace_back(x, y, m.second, start_score);

        // Each entry in the tree points to an index in the node list
        tree.insert(x, y, nodes.size() - 1);
    }

    // Construct graph based on nearest down-diagonal neighbors in the spatial index
    for (size_t i=0; i<nodes.size(); i++){
        if (i <= SINK_INDEX){
            continue;
        }

        auto& node = nodes[i];

        size_t x_stop = node.get_x_stop();
        size_t y_stop = node.get_y_stop();

        vector <pair <coord_t, size_t> > results;
        tree.find(x_stop, y_stop, max_gap, results, true);

        // Given the nearest down-diagonal matches within a certain gap size, add edges in the graph to those matches
        for (const auto& r: results){
            auto& adjacent_node = nodes[r.second];
            node.next.emplace_back(r.second);
            adjacent_node.prev.emplace_back(i);
            n_edges++;

//            cerr << node.x << ' ' << node.y << " (" << r.first.first << ',' << r.first.second << ")=" << r.second << '\n';
        }

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


size_t Dag::get_x_diagonal(Node& n){
    return (y_size - 1) - n.y + n.x;
}


size_t Dag::get_y_diagonal(Node& n){
    return n.x + n.y;
}


void Dag::step_alignment_forward(size_t i, vector<size_t>& max_prev){
    if (i == SOURCE_INDEX){
        nodes[i].score = 0;
    }
    else {
        for (const auto j: nodes[i].prev){
            int64_t prev_diagonal = get_x_diagonal(nodes[j]);
            int64_t diagonal = get_x_diagonal(nodes[i]);
            int64_t gap_size = abs(diagonal-prev_diagonal);

            int64_t score = nodes[j].score + nodes[i].length*match_score;

            // No gap cost for source/sink gaps, and no gap_open if gap=0
            if (gap_size > 0 and j != SOURCE_INDEX and i != SINK_INDEX){
                score += gap_size*gap_score + gap_open_score;
            }

            if (score > nodes[i].score){
                nodes[i].score = score;
                max_prev[i] = j;
            }
        }
    }
}


void Dag::compute_alignment(){
    vector <size_t> max_prev(nodes.size(),0);

    // Forward pass
    for_each_in_kahns_iteration([&](const size_t i){
        step_alignment_forward(i, max_prev);
    });

    // Traceback
    auto i = SINK_INDEX;
    while (i != SOURCE_INDEX){
        alignment_path.emplace_front(i);
        i = max_prev[i];
    }
}


void Dag::write_to_svg(path output_path, bool autoscale){
    cerr << "Writing SVG: " << output_path << '\n';

    size_t width_px = 800;
    size_t height_px = 800;

    size_t line_width;
    if (autoscale) {
        line_width = max(width_px,height_px)*0.5;
    }
    else{
        line_width = 30;
    }

    SvgPlot plot(output_path, width_px, height_px, 0, x_size, 0, y_size, true);
    string point_type = "circle";

    // Plot all the matches and their edges
    for (size_t i=0; i<nodes.size(); i++){
        auto& node = nodes[i];

        size_t x_stop = node.get_x_stop();
        size_t y_stop = node.get_y_stop();

        // Draw the match element
        plot.add_line(node.x, node.y, x_stop, y_stop, double(line_width)/120.0, "gray");
        plot.add_point(node.x, node.y, point_type, double(line_width)/120.0, "black");

        for (size_t j: node.prev){
            auto& prev_node = nodes[j];

            auto prev_x_stop = prev_node.get_x_stop();
            auto prev_y_stop = prev_node.get_y_stop();

            // Draw edge between matches
            plot.add_curve(prev_x_stop, prev_y_stop, node.x, node.y, double(line_width)/40.0, double(line_width)/800.0, "orange");
        }
    }

    // Plot the alignment path that was computed if it exists
    for (size_t path_index=0; path_index<alignment_path.size(); path_index++){
        auto i = alignment_path[path_index];

        auto& node = nodes[i];

        size_t x_stop = node.get_x_stop();
        size_t y_stop = node.get_y_stop();

        // Overlay on top of match
        plot.add_line(node.x, node.y, x_stop, y_stop, double(line_width)/600.0, "red");
        plot.add_point(node.x, node.y, point_type, double(line_width)/180.0, "red");

        if (path_index > 0){
            auto j = alignment_path[path_index-1];

            auto& prev_node = nodes[j];

            auto prev_x_stop = prev_node.get_x_stop();
            auto prev_y_stop = prev_node.get_y_stop();

            // Overlay on top of edge
            plot.add_curve(prev_x_stop, prev_y_stop, node.x, node.y, double(line_width)/40.0, double(line_width)/1000.0, "red");
            plot.add_point(node.x, node.y, point_type, double(line_width)/180.0, "red");
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
