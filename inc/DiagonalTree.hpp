#ifndef OVERLAP_ANALYSIS_DIAGONALTREE_HPP
#define OVERLAP_ANALYSIS_DIAGONALTREE_HPP

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



class DiagonalTree {
public:
    vector <map <size_t, size_t> > diagonals;
    size_t x_size;
    size_t y_size;

    DiagonalTree(size_t x_size, size_t y_size);
    DiagonalTree(const vector <pair <pair <size_t,size_t>, size_t> >& matches, size_t x_size, size_t y_size);

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


DiagonalTree::DiagonalTree(size_t x_size, size_t y_size):
    diagonals(x_size + y_size - 1),
    x_size(x_size),
    y_size(y_size)
{}


DiagonalTree::DiagonalTree(const vector <pair <pair <size_t,size_t>, size_t> >& matches, size_t x_size, size_t y_size):
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


size_t DiagonalTree::get_x_diagonal(size_t x, size_t y){
//    if (x >= x_size or y >= y_size){
//        throw runtime_error("ERROR: cannot compute diagonal for " + to_string(x) + ',' + to_string(y) +
//                            " because >= matrix size: " + to_string(x_size) + ',' + to_string(y_size));
//    }
    return (y_size - 1) - y + x;
}


size_t DiagonalTree::get_y_diagonal(size_t x, size_t y){
    return x + y;
}


size_t DiagonalTree::get_x(size_t x_diag, size_t y_diag){
    return (y_diag + x_diag + 1 - y_size) / 2;
}


size_t DiagonalTree::get_y(size_t x_diag, size_t y_diag){
    return (y_diag - x_diag + y_size - 1) / 2;
}


void DiagonalTree::insert(size_t x, size_t y, size_t value){
    size_t x_diag = get_x_diagonal(x,y);
    size_t y_diag = get_y_diagonal(x,y);

    diagonals[x_diag].emplace(y_diag, value);
}


void DiagonalTree::get_range_from_tree(
        size_t x_diag,
        size_t a,
        size_t b,
        vector<pair <pair <size_t, size_t>, size_t> >& results,
        bool first_only){

    auto lower = diagonals[x_diag].lower_bound(a);
    auto upper = diagonals[x_diag].upper_bound(b);

    if (lower != upper){
        for (auto it=lower; it!=upper; ++it){
            auto x = get_x(x_diag, it->first);
            auto y = get_y(x_diag, it->first);

            results.emplace_back(make_pair(x,y), it->second);

            if (first_only){
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
void DiagonalTree::find(size_t x, size_t y, size_t box_length, vector<pair <pair <size_t, size_t>, size_t> >& results, bool first_only) {
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

#endif //OVERLAP_ANALYSIS_DIAGONALTREE_HPP
