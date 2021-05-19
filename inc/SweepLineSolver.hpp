#ifndef OVERLAP_ANALYSIS_SWEEPLINESOLVER_HPP
#define OVERLAP_ANALYSIS_SWEEPLINESOLVER_HPP

#include <string>
#include <vector>
#include <utility>

using std::string;
using std::vector;
using std::pair;


class Terminus {
public:
    size_t diagonal;
    size_t antidiagonal;
    size_t index;
    bool is_start;
    bool is_removal_flag;

    Terminus(size_t diagonal, size_t antidiagonal, size_t match_index, bool is_start, bool is_removal_flag);
};


class DiagonalMatch {
public:
    /// Attributes ///
    size_t diagonal;
    size_t antidiagonal;
    size_t length;

    DiagonalMatch* nearest_adjacent;
    DiagonalMatch* nearest_down_diagonal;

    int64_t adjacent_gap_size = -1;
    int64_t down_diagonal_gap_size = -1;

    /// Methods ///
    DiagonalMatch(size_t diagonal, size_t antidiagonal, size_t length);
    DiagonalMatch(size_t x, size_t y, size_t length, size_t y_size);

    size_t get_x_diagonal(size_t x, size_t y, size_t y_size) const;
    size_t get_y_diagonal(size_t x, size_t y) const;
    size_t get_x(size_t y_size) const;
    size_t get_y(size_t y_size) const;
};


class SweepLineSolver {
public:
    /// Attributes ///
    size_t diagonal_size;
    size_t antidiagonal_size;

    size_t x_size;
    size_t y_size;

    int64_t max_gap = 1;
    int64_t max_skip = 3;

    /// Methods ///
    SweepLineSolver(vector <DiagonalMatch>& matches);

    bool is_adjacent(const Terminus& t1, const Terminus& t2);

    void try_add_edge(vector <DiagonalMatch>& matches,
                      const Terminus& t1,
                      const Terminus& t2,
                      int64_t& gap_size,
                      int64_t& skip_size);
};




#endif //OVERLAP_ANALYSIS_SWEEPLINESOLVER_HPP
