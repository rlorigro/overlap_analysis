#ifndef OVERLAP_ANALYSIS_DIAGONALINTERVALTREE_HPP
#define OVERLAP_ANALYSIS_DIAGONALINTERVALTREE_HPP

#include <stdexcept>
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <tuple>
#include <cmath>

using std::runtime_error;
using std::to_string;
using std::cerr;
using std::map;
using std::vector;
using std::pair;
using std::make_pair;
using std::tie;
using std::min;

#include "boost/icl/interval_set.hpp"
#include "boost/icl/interval.hpp"

using boost::icl::interval_set;
using boost::icl::interval;


class DiagonalIntervalTree {
public:
    vector <interval_set <uint32_t> > diagonals;
    size_t x_size;
    size_t y_size;

    DiagonalIntervalTree(size_t x_size, size_t y_size);
    DiagonalIntervalTree(const vector <pair <pair <size_t,size_t>, size_t> >& matches, size_t x_size, size_t y_size);

    size_t get_x_diagonal(size_t x, size_t y);
    size_t get_y_diagonal(size_t x, size_t y);

    size_t get_x(size_t x_diag, size_t y_diag);
    size_t get_y(size_t x_diag, size_t y_diag);

    void insert(size_t x, size_t y, size_t value);

    void find(
            size_t x,
            size_t y,
            size_t box_length,
            vector<pair <pair <size_t, size_t>, size_t> >& results,
            bool first_only=false);

    void get_range_from_tree(
            size_t x_diag,
            size_t a,
            size_t b,
            vector<pair <pair <size_t, size_t>, size_t> >& results,
            bool first_only=false);
};


DiagonalIntervalTree::DiagonalIntervalTree(size_t x_size, size_t y_size):
    diagonals(x_size + y_size - 1),
    x_size(x_size),
    y_size(y_size)
{}


DiagonalIntervalTree::DiagonalIntervalTree(const vector <pair <pair <size_t,size_t>, size_t> >& matches, size_t x_size, size_t y_size):
    diagonals(x_size + y_size - 1),
    x_size(x_size),
    y_size(y_size)
{
    for (const auto& m: matches){
        auto x = m.first.first;
        auto y = m.first.second;

        insert(x, y, m.second);
    }
}


size_t DiagonalIntervalTree::get_x_diagonal(size_t x, size_t y){
//    if (x >= x_size or y >= y_size){
//        throw runtime_error("ERROR: cannot compute diagonal for " + to_string(x) + ',' + to_string(y) +
//                            " because >= matrix size: " + to_string(x_size) + ',' + to_string(y_size));
//    }
    return (y_size - 1) - y + x;
}


size_t DiagonalIntervalTree::get_y_diagonal(size_t x, size_t y){
    return x + y;
}


size_t DiagonalIntervalTree::get_x(size_t x_diag, size_t y_diag){
    return (y_diag + x_diag + 1 - y_size) / 2;
}


size_t DiagonalIntervalTree::get_y(size_t x_diag, size_t y_diag){
    return (y_diag - x_diag + y_size - 1) / 2;
}


void DiagonalIntervalTree::insert(size_t x, size_t y, size_t length){
    size_t x_diag = get_x_diagonal(x,y);
    size_t y_diag = get_y_diagonal(x,y);

    // Construct a Boost interval for this numeric range/coord
    auto a = interval<uint32_t>::closed(y_diag, y_diag+length-1);

    diagonals[x_diag].add(a);
}


void DiagonalIntervalTree::get_range_from_tree(
        size_t x_diag,
        size_t a,
        size_t b,
        vector<pair <pair <size_t, size_t>, size_t> >& results,
        bool first_only){

    // Construct a Boost interval for this query range/coord
    auto range = interval<uint32_t>::closed(a, b);

    auto lower = diagonals[x_diag].lower_bound(range);
    auto upper = diagonals[x_diag].upper_bound(range);

    if (lower != upper){
        for (auto it=lower; it!=upper; ++it){
            // Convert the interval result start coordinate back into x,y space
            auto x = get_x(x_diag, it->lower());
            auto y = get_y(x_diag, it->lower());

            auto length = it->upper() - it->lower() + 1;

            results.emplace_back(make_pair(x,y), length);

            // Stop after the first result on this diagonal that starts within the query region [a,b]
            if (first_only and it->lower() >= a){
                break;
            }
        }
    }
}


///
/// \param x
/// \param y
/// \param box_length
/// \param results
/// \param first_only If true, then only the first result per diagonal will be returned.
///
//     Box length specifies the region beyond the given xy point to search. For example, if box length = 1:
//
//         x 1
//        ._._.
//      y |_|_|
//      1 |_|_|
//
//      This entire grid will return results if any are contained.
//
void DiagonalIntervalTree::find(size_t x, size_t y, size_t box_length, vector<pair <pair <size_t, size_t>, size_t> >& results, bool first_only) {
    if (x > x_size or y > y_size){
        return;
    }

    size_t search_space;

    search_space = box_length + 1;

    // Iterate the top edge of the box assuming matrix orientation s.t. (0,0) is top left
    for (size_t i=x; i<=x+box_length; i++){
        if (i >= x_size){
            continue;
        }

        auto x_diag = get_x_diagonal(i,y);
        auto y_diag = get_y_diagonal(i,y);

        get_range_from_tree(x_diag, y_diag, y_diag+search_space, results, first_only);
        search_space--;
    }

    search_space = box_length;

    // Iterate the left edge of the box assuming matrix orientation s.t. (0,0) is top left
    // Don't re-iterate the corner of the box so start at y+1
    for (size_t i=y+1; i<=y+box_length; i++){
        if (i >= y_size){
            continue;
        }

        auto x_diag = get_x_diagonal(x,i);
        auto y_diag = get_y_diagonal(x,i);

        get_range_from_tree(x_diag, y_diag, y_diag+search_space, results, first_only);
        search_space--;
    }
}

#endif //OVERLAP_ANALYSIS_DIAGONALINTERVALTREE_HPP
