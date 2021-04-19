#include "DagAligner.hpp"

using std::make_shared;


Node::Node(size_t x, size_t y, size_t length):
    x(x),
    y(y),
    length(length)
{}


Dag::Dag(const vector<xy_match>& matches, size_t max_gap):
        tree(matches),
        nodes(),
        source(make_shared<Node>(Node(0,0,0))),
        sink(make_shared<Node>(Node(0,0,0))),
        max_gap(max_gap)
{
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
            cerr << x_start << ' ' << y_start << " (" << r.first.get<0>() << ',' << r.first.get<1>() << ")=" << r.second << '\n';
        }
        cerr << '\n';
    }
}