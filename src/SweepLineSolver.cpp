#include "SweepLineSolver.hpp"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <set>
#include <map>

using std::runtime_error;
using std::sort;
using std::to_string;
using std::cerr;
using std::set;
using std::map;


DiagonalMatch::DiagonalMatch(size_t diagonal, size_t antidiagonal, size_t length):
    diagonal(diagonal),
    antidiagonal(antidiagonal),
    length(length),
    nearest_adjacent(nullptr),
    nearest_down_diagonal(nullptr)
{}


DiagonalMatch::DiagonalMatch(size_t x, size_t y, size_t length, size_t y_size):
    diagonal(get_x_diagonal(x, y, y_size)),
    antidiagonal(get_y_diagonal(x, y)),
    length(length),
    nearest_adjacent(nullptr),
    nearest_down_diagonal(nullptr)
{}


size_t DiagonalMatch::get_x_diagonal(size_t x, size_t y, size_t y_size) const{
    return (y_size - 1) - y + x;
}


size_t DiagonalMatch::get_y_diagonal(size_t x, size_t y) const{
    return x + y;
}


size_t DiagonalMatch::get_x(size_t y_size) const{
    return (antidiagonal + diagonal + 1 - y_size) / 2;
}


size_t DiagonalMatch::get_y(size_t y_size) const{
    return (antidiagonal - diagonal + y_size - 1) / 2;
}


bool operator<(DiagonalMatch& a, DiagonalMatch& b){
    return a.antidiagonal < b.antidiagonal;
}


bool operator<(Terminus& a, Terminus& b){
    return a.diagonal < b.diagonal;
}


bool operator<(const Terminus& a, const Terminus& b){
    return a.diagonal < b.diagonal;
}


Terminus::Terminus(size_t diagonal, size_t antidiagonal, size_t index, bool is_start):
    diagonal(diagonal),
    antidiagonal(antidiagonal),
    index(index),
    is_start(is_start)
{}


SweepLineSolver::SweepLineSolver(const vector <DiagonalMatch>& matches){
    vector<Terminus> termini;

    // Sort all the start and end positions of the matches (line segments) in the vector.
    for (size_t i=0; i<matches.size(); i++){
        const auto& m = matches[i];
        termini.emplace_back(m.diagonal, m.antidiagonal, i, true);
        termini.emplace_back(m.diagonal, m.antidiagonal + 2*m.length, i, false);
    }

    sort(termini.begin(), termini.end(), [&](const Terminus& a, const Terminus& b){
        if (a.antidiagonal < b.antidiagonal){
            return true;
        }
        // If 2 termini are in the same coordinate, then the one that is ending goes first
        else if (a.antidiagonal == b.antidiagonal){
            return !a.is_start;
        }
        else{
            return false;
        }
    });

    int64_t cutoff = 1;

    set <Terminus> visited;
    for (auto& item: termini){
        if (item.is_start) {
            // Do something
        }
        else{
            // Do something else
        }

        cerr << "inserting " << item.diagonal << '\n';

        auto result = visited.emplace(item);
        bool success = result.second;
        auto iter = result.first;

        if (not success){
            // If this is the end terminus of a previously existing start terminus, do something special
            if (iter->index == item.index and iter->is_start and not item.is_start){
                cerr << "ending match" << '\n';
            }
            else {
                cerr << "ERROR: element in sweep line has conflict with existing element on diagonal" << '\n';
            }
        }
        else{
            // Anything need to be here?
        }

        size_t n_steps;
        int64_t distance;

        n_steps = 0;
        distance = 0;
        cerr << "Checking lower bounds" << '\n';
        while (++iter != visited.end() and distance <= cutoff){
            distance = llabs(int64_t(iter->diagonal) - int64_t(item.diagonal));
            n_steps++;

            if (distance <= cutoff){
                cerr  << item.diagonal << " " << iter->diagonal << '\n';
            }
        }

        std::advance(iter,-n_steps);

        n_steps = 0;
        distance = 0;
        cerr << "Checking upper bounds" << '\n';
        while (iter-- != visited.begin() and distance <= cutoff){
            distance = llabs(int64_t(iter->diagonal) - int64_t(item.diagonal));
            n_steps++;

            if (distance <= cutoff){
                cerr  << item.diagonal << " " << iter->diagonal << '\n';
            }
        }

        cerr << '\n';
    }
}


