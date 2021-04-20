#ifndef OVERLAP_ANALYSIS_DAGALIGNER_H
#define OVERLAP_ANALYSIS_DAGALIGNER_H

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

using boost::geometry::wkt;
using boost::geometry::model::point;
using boost::geometry::model::box;
using boost::geometry::cs::cartesian;
using boost::geometry::index::rtree;
using boost::geometry::index::quadratic;
using boost::geometry::index::intersects;

#include <experimental/filesystem>
#include <stdexcept>
#include <cmath>
#include <chrono>
#include <vector>
#include <iostream>
#include <memory>
#include <utility>
#include <queue>
#include <deque>
#include <iostream>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::path;
using std::runtime_error;
using std::shared_ptr;
using std::vector;
using std::cerr;
using std::make_pair;
using std::pair;
using std::queue;
using std::deque;
using std::ostream;


using coord_t = pair<size_t,size_t>;


// Required data type for using the Boost spatial indexing data structure
using xy_point = point<size_t, 2, cartesian>;

// Represents an exact match found in a 2d alignment matrix, with x=ref, y=query, and its length along the diagonal
using xy_match = pair<xy_point, size_t>;

// Graph form of an exact match found in a 2d alignment matrix, with x=ref, y=query, and its length along the diagonal.
// x,y coords are used to compute gap (edge) scores later, based on some scoring scheme
class Node {
public:
    // vector indexes which point to neighbor nodes
     vector <size_t> prev;
     vector <size_t> next;

    size_t x;
    size_t y;
    size_t length;
    int64_t score;

    Node(size_t x, size_t y, size_t length);
};


ostream& operator<<(ostream& o, Node& n);


class Dag {
public:
    vector <Node> nodes;
    deque <size_t> alignment_path;
    size_t x_size;
    size_t y_size;
    size_t n_edges;
    size_t max_gap;

    const size_t SOURCE_INDEX = 0;
    const size_t SINK_INDEX = 1;
    const size_t START_INDEX = 2;

    Dag(const vector<xy_match>& matches, size_t x_size, size_t y_size, size_t max_gap);
    Dag(const vector <pair <coord_t,size_t> >& matches, size_t x_size, size_t y_size, size_t max_gap);

    template <class T> void for_each_in_kahns_iteration(T function);
    void compute_alignment();

    void write_to_svg(path output_path);
};


//  Kahns algorithm (adapted to be non-destructive from https://en.wikipedia.org/wiki/Topological_sorting)
//  ------------------------------------------------------------------------------------------------------
//      L ← Empty list that will contain the sorted elements
//      S ← Set of all nodes with no incoming edge
//
//      while S is not empty do
//        remove a node n from S
//        add n to L
//        for each node m with an edge e from n to m do
//            remove edge e from the graph
//            if m has no other incoming edges then
//                  insert m into S
//
//      if graph has edges then
//        return error   (graph has at least one cycle)
//      else
//        return L   (a topologically sorted order)
//
template <class T> void Dag::for_each_in_kahns_iteration(T f){
    vector<size_t> edge_counts(nodes.size(),0);
    size_t edges_counted = 0;
    vector<size_t> l;
    queue<size_t> s;

    s.push(SOURCE_INDEX);

    while (not s.empty()){
        const auto i = s.front();
        s.pop();

        // Where the work is done, at this step it is guaranteed
        // that all previous paths have been visited
        f(i);
        l.emplace_back(i);

        cerr << i << ' ' << nodes[i].next.size() << '\n';

        for (auto j: nodes[i].next){
            cerr << "next node: " << j << '\n';

            edge_counts[j]++;
            edges_counted++;

            if (nodes[j].prev.size() == edge_counts[j]){
                s.push(j);
            }
        }
    }

    if (edges_counted != n_edges){
        throw runtime_error("ERROR: graph contains cycles");
    }
}

ostream& operator<<(ostream& o, Dag& dag);

#endif //OVERLAP_ANALYSIS_DAGALIGNER_H
