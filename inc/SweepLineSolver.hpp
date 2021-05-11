#ifndef OVERLAP_ANALYSIS_SWEEPLINESOLVER_HPP
#define OVERLAP_ANALYSIS_SWEEPLINESOLVER_HPP

#include <vector>
#include <utility>

using std::vector;
using std::pair;


class Terminus {
public:
    size_t diagonal;
    size_t antidiagonal;
    size_t index;
    bool is_start;

    Terminus(size_t diagonal, size_t antidiagonal, size_t match_index, bool is_start);
};


class DiagonalMatch {
public:
    /// Attributes ///
    size_t diagonal;
    size_t antidiagonal;
    size_t length;

    DiagonalMatch* nearest_adjacent;
    DiagonalMatch* nearest_down_diagonal;

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

    /// Methods ///
    SweepLineSolver(const vector <DiagonalMatch>& matches);
};




#endif //OVERLAP_ANALYSIS_SWEEPLINESOLVER_HPP
