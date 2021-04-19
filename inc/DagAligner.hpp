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

#include <cmath>
#include <chrono>
#include <vector>
#include <iostream>
#include <memory>
#include <utility>

using std::shared_ptr;
using std::vector;
using std::cerr;
using std::make_pair;
using std::pair;


// Required data type for using the Boost spatial indexing data structure
using xy_point = point<size_t, 2, cartesian>;

// Represents an exact match found in a 2d alignment matrix, with x=ref, y=query, and its length along the diagonal
using xy_match = pair<xy_point, size_t>;

// Graph form of an exact match found in a 2d alignment matrix, with x=ref, y=query, and its length along the diagonal.
// x,y coords are used to compute gap (edge) scores later, based on some scoring scheme
class Node {
public:
    vector <shared_ptr <Node> > prev;
    vector <shared_ptr <Node> > next;
    size_t x;
    size_t y;
    size_t length;

    Node(size_t x, size_t y, size_t score);
};


class Dag {
public:
    rtree <xy_match, quadratic<30> > tree;
    vector <shared_ptr <Node> > nodes;
    shared_ptr<Node> source;
    shared_ptr<Node> sink;
    size_t max_gap;

    Dag(const vector<xy_match>& matches, size_t max_gap);
//    void add_edge(shared_ptr<Node> a, shared_ptr<Node> b);
};




#endif //OVERLAP_ANALYSIS_DAGALIGNER_H
