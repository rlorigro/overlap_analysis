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


Terminus::Terminus(size_t diagonal, size_t antidiagonal, size_t index, bool is_start, bool is_removal_flag):
    diagonal(diagonal),
    antidiagonal(antidiagonal),
    index(index),
    is_start(is_start),
    is_removal_flag(is_removal_flag)
{}


bool SweepLineSolver::is_adjacent(const Terminus& t1, const Terminus& t2){
    return t1.is_start and t2.is_start;
}


void SweepLineSolver::try_add_edge(vector <DiagonalMatch>& matches, const Terminus& t1, const Terminus& t2, int64_t& gap_size, int64_t& skip_size){
    gap_size = llabs(int64_t(t1.diagonal) - int64_t(t2.diagonal));
    skip_size = llabs(int64_t(t1.antidiagonal) - int64_t(t2.antidiagonal));

    if (gap_size <= max_gap and skip_size < max_skip){
        cerr << "  Creating edge: (" << t2.diagonal << ',' << t2.antidiagonal << ")->(" << t1.diagonal << ',' << t1.antidiagonal << ")" << '\n';

        auto& match1 = matches[t1.index];
        auto& match2 = matches[t2.index];

        bool adjacent = is_adjacent(t1,t2);

        // Matches overlap
        if (adjacent){
            if (gap_size < match1.adjacent_gap_size){
                match1.nearest_adjacent = &match2;
            }
            if (gap_size < match2.adjacent_gap_size){
                match2.nearest_adjacent = &match1;
            }
        }
        // Match2 is up-diagonal from match1
        else if (t1.is_start and not t2.is_start){
            if (gap_size < match2.down_diagonal_gap_size or match2.nearest_down_diagonal == nullptr){
                match2.nearest_down_diagonal = &match1;
            }
        }
        // Match2 is down-diagonal from match1
        else if (t2.is_start and not t1.is_start){
            if (gap_size < match1.down_diagonal_gap_size or match1.nearest_down_diagonal == nullptr){
                match1.nearest_down_diagonal = &match2;
            }
        }
        // Anything else is not a valid condition for creating an edge
    }
}


SweepLineSolver::SweepLineSolver(vector <DiagonalMatch>& matches){
    vector<Terminus> termini;

    for (size_t i=0; i<matches.size(); i++){
        const auto& m = matches[i];
        termini.emplace_back(m.diagonal, m.antidiagonal, i, true, false);
        termini.emplace_back(m.diagonal, m.antidiagonal + 2*m.length, i, false, false);
        termini.emplace_back(m.diagonal, m.antidiagonal + 2*m.length + max_skip, i, false, true);
    }

    // Sort all the start, end, and end-of-scope positions of the matches (line segments) in the vector.
    sort(termini.begin(), termini.end(), [&](const Terminus& a, const Terminus& b){
        if (a.antidiagonal < b.antidiagonal){
            return true;
        }
        // If 2 termini are in the same coordinate, then the one that is ending goes first
        else if (a.antidiagonal == b.antidiagonal){
            int32_t a_ordinal = !a.is_start + a.is_removal_flag;
            int32_t b_ordinal = !b.is_start + b.is_removal_flag;

            return a_ordinal < b_ordinal;
        }
        else{
            return false;
        }
    });

    set <Terminus> visited;
    for (auto& terminus: termini){
        if (terminus.is_start) {
            // Do something
        }
        else{
            // Do something else
        }

        cerr << '\n';
        cerr << "--- inserting " << terminus.diagonal << '\n';

        bool success;
        set<Terminus>::iterator other_terminus;
        tie(other_terminus, success) = visited.emplace(terminus);

        int64_t gap_size = 0;
        int64_t skip_size = 0;
        size_t n_steps;

        if (not success){
            // Check if this entry has a conflict with a previously existing entry
            if (other_terminus->index == terminus.index){
                // If this is the ending terminus, simply replace the element
                if (other_terminus->is_start and not terminus.is_start){
                    cerr << "--- ending match" << '\n';

                    visited.erase(other_terminus);
                    tie(other_terminus, success) = visited.emplace(terminus);

                    if (not success){
                        throw runtime_error("ERROR: conflicting elements cannot be resolved");
                    }
                }
                // If this terminus is a flag indicating the end-of-scope of a previous terminus, simply erase previous
                else if (not other_terminus->is_start and terminus.is_removal_flag){
                    cerr << "--- erasing out-of-scope terminus" << '\n';
                    visited.erase(other_terminus);
                }
                else {
                    throw runtime_error("ERROR: new element has same match index as terminated element");
                }
            }
            else {
                // If the end terminus is being overwritten by a new terminus, attempt to add an edge before overwriting
                if (not other_terminus->is_start and terminus.is_start){
                    cerr << "--- overwriting end terminus with start terminus" << '\n';
                    try_add_edge(matches, *other_terminus, terminus, gap_size, skip_size);

                    visited.erase(other_terminus);
                    tie(other_terminus, success) = visited.emplace(terminus);

                    if (not success){
                        throw runtime_error("ERROR: conflicting elements cannot be resolved");
                    }
                }
                // Ignore removal flag if the terminus has already been replaced by something
                else if (terminus.is_removal_flag){
                    // Do nothing
                }
                else{
                    throw runtime_error("ERROR: new element conflicts with unterminated previous element");
                }
            }
        }
        else{
            // Anything need to be here?
        }

        if (terminus.is_removal_flag){
            continue;
        }

        n_steps = 0;
        gap_size = 0;
        cerr << "Checking lower bounds" << '\n';
        while (++other_terminus != visited.end() and gap_size <= max_gap){
            try_add_edge(matches, *other_terminus, terminus, gap_size, skip_size);
            n_steps++;
        }

        // Skip back to the other side of the tree
        std::advance(other_terminus, -(n_steps+1));

        n_steps = 0;
        gap_size = 0;
        cerr << "Checking upper bounds" << '\n';
        while (other_terminus-- != visited.begin() and gap_size <= max_gap){
            try_add_edge(matches, *other_terminus, terminus, gap_size, skip_size);
            n_steps++;
        }
    }
}


