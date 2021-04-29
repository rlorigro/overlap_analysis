//
//  wavefront.cpp
//  
//
//  Created by Jordan Eizenga on 4/4/21.
//

#include "wavefront.hpp"
#include "SvgPlot.hpp"
#include "Color.hpp"

#include <algorithm>
#include <iostream>
#include <deque>

#include <cmath>

using std::round;
using std::min;
using std::max;

struct Wavefront {
public:
    Wavefront() {}
    Wavefront(int32_t diag_begin, int32_t diag_end)
            : diag_begin(diag_begin), M(diag_end - diag_begin,
                                        std::numeric_limits<int32_t>::min()), D(M), I(M) {

    }

    int32_t diag_begin;
    std::deque<int32_t> M, I, D;
};

void wavefront_extend(const std::string& seq1, const std::string& seq2,
                      Wavefront& wf) {

    for (int64_t k = 0; k < wf.M.size(); ++k) {
        int64_t diag = wf.diag_begin + k;
        int64_t anti_diag = wf.M[k];
        int64_t i = (diag + anti_diag) / 2 + 1;
        int64_t j = (anti_diag - diag) / 2 + 1;
        while (i < seq1.size() && j < seq2.size() && seq1[i] == seq2[j]) {
            wf.M[k] += 2;
            ++i;
            ++j;
        }
    }
}

void wavefront_next(const std::string& seq1, const std::string& seq2,
                    const WFScores& scores, std::vector<Wavefront>& wfs) {


    // find the range of incoming wavefronts
    int32_t hi = std::numeric_limits<int32_t>::min();
    int32_t lo = std::numeric_limits<int32_t>::max();
    for (auto s_back : {scores.mismatch, scores.gap_extend, scores.gap_open + scores.gap_extend}) {
        if (wfs.size() >= s_back) {
            const auto& wf_prev = wfs[wfs.size() - s_back];
            if (!wf_prev.M.empty()) {
                lo = std::min<int32_t>(lo, wf_prev.diag_begin);
                hi = std::max<int32_t>(hi, wf_prev.diag_begin + (int64_t) wf_prev.M.size());
            }
        }
    }

    if (hi < lo) {
        // there are no incoming wavefronts
        wfs.emplace_back(0, 0);
        return;
    }

    // expand the range by 1
    hi += 1;
    lo -= 1;

    wfs.emplace_back(lo, hi);
    auto& wf = wfs.back();

    // open gaps
    if (wfs.size() >= scores.gap_extend + scores.gap_open + 1) {
        const auto& wf_prev = wfs[wfs.size() - scores.gap_extend - scores.gap_open - 1];
        for (auto k = lo; k < hi; ++k) {
            if (k - 1 >= wf_prev.diag_begin && k - 1 < wf_prev.diag_begin + (int64_t) wf_prev.M.size()) {
                int64_t a = wf_prev.M[k - 1 - wf_prev.diag_begin] + 1;
                if ((k + a) / 2 < (int64_t) seq1.size() && (a - k) / 2 < (int64_t) seq2.size()) {
                    wf.I[k - lo] = a;
                }
            }
            if (k + 1 >= wf_prev.diag_begin && k + 1 < wf_prev.diag_begin + (int64_t) wf_prev.M.size()) {
                int64_t a = wf_prev.M[k + 1 - wf_prev.diag_begin] + 1;
                if ((k + a) / 2 < (int64_t) seq1.size() && (a - k) / 2 < (int64_t) seq2.size()) {
                    wf.D[k - lo] = a;
                }
            }
        }
    }

    // extend gaps
    if (wfs.size() >= scores.gap_extend + 1) {
        const auto& wf_prev = wfs[wfs.size() - scores.gap_extend - 1];
        for (auto k = lo; k < hi; ++k) {
            if (k - 1 >= wf_prev.diag_begin && k - 1 < wf_prev.diag_begin + (int64_t) wf_prev.M.size()) {
                int64_t a = wf_prev.I[k - 1 - wf_prev.diag_begin] + 1;
                if ((k + a) / 2 < (int64_t) seq1.size() && (a - k) / 2 < (int64_t) seq2.size()) {
                    wf.I[k - lo] = std::max<int32_t>(wf.I[k - lo], a);
                }
            }
            if (k + 1 >= wf_prev.diag_begin && k + 1 < wf_prev.diag_begin + (int64_t) wf_prev.M.size()) {
                int64_t a = wf_prev.D[k + 1 - wf_prev.diag_begin] + 1;
                if ((k + a) / 2 < (int64_t) seq1.size() && (a - k) / 2 < (int64_t) seq2.size()) {
                    wf.D[k - lo] = std::max<int32_t>(wf.D[k - lo], a);
                }
            }
        }
    }


    // mismatches
    if (wfs.size() >= scores.mismatch + 1) {
        const auto& wf_prev = wfs[wfs.size() - scores.mismatch - 1];
        for (auto k = lo; k < hi; ++k) {
            if (k >= wf_prev.diag_begin && k < wf_prev.diag_begin + (int64_t) wf_prev.M.size()) {
                int64_t a = wf_prev.M[k - wf_prev.diag_begin] + 2;
                if ((k + a) / 2 < (int64_t) seq1.size() && (a - k) / 2 < (int64_t) seq2.size()) {
                    wf.M[k - lo] = a;
                }
            }
        }
    }

    // calculate the opt
    for (size_t k = 0; k < wf.M.size(); ++k) {
        wf.M[k] = std::max(wf.M[k], std::max(wf.I[k], wf.D[k]));
    }
}

bool wavefront_reached(const Wavefront& wf, int32_t diag, int32_t anti_diag) {
    for (auto k = 0; k < wf.M.size(); ++k) {
        if (wf.diag_begin + k == diag && wf.M[k] == anti_diag) {
            return true;
        }
    }
    return false;
}

void wavefront_traceback(const std::string& seq1, const std::string& seq2,
                         const WFScores& scores, const std::vector<Wavefront>& wfs,
                         std::vector<CIGAROp>& cigar) {

    enum WFMatrix_t {MAT_M, MAT_D, MAT_I};

    int64_t d = seq1.size() - seq2.size();
    int64_t a = seq1.size() + seq2.size() - 2;
    int64_t s = wfs.size() - 1;
    WFMatrix_t mat = MAT_M;

    uint32_t op_len = 0;
    while (a > -2) {
        // we're in the match/mismatch matrix
        const auto& wf = wfs[s];
        if (mat == MAT_M) {
            // extend through any matches
            int64_t i = (d + a) / 2;
            int64_t j = (a - d) / 2;
            while (i >= 0 && j >= 0 && seq1[i] == seq2[j]
                   && a != wf.I[d - wf.diag_begin] && a != wf.D[d - wf.diag_begin]) {
                op_len += 1;
                --i;
                --j;
                a -= 2;
            }
            if (a == -2) {
                // this is the beginning of the alignment, so don't
                // look for the score-increasing edge before this run
                break;
            }
            if (a == wf.I[d - wf.diag_begin]) {
                // this is where an insertion closed
                cigar.emplace_back('M', op_len);
                mat = MAT_I;
                op_len = 0;
                continue;
            }
            else if (a == wf.D[d - wf.diag_begin]) {
                // this is where a deletion closed
                cigar.emplace_back('M', op_len);
                mat = MAT_D;
                op_len = 0;
                continue;
            }
            else {
                // this is a mismatch
                op_len += 1;
                a -= 2;
                s -= scores.mismatch;
                continue;
            }
        }
        else if (mat == MAT_I) {
            // we're in the insertion matrix
            if (s >= scores.gap_extend) {
                const auto& wf_prev = wfs[s - scores.gap_extend];
                if (d - 1 >= wf_prev.diag_begin && d - 1 < wf_prev.diag_begin + (int64_t) wf_prev.M.size()) {
                    if (wf.I[d - wf.diag_begin] == wf_prev.I[d - 1 - wf_prev.diag_begin] + 1) {
                        // an insert extended here
                        s -= scores.gap_extend;
                        op_len += 1;
                        d -= 1;
                        a -= 1;
                        continue;
                    }
                }
            }
            if (s >= scores.gap_extend + scores.gap_open) {
                const auto& wf_prev = wfs[s - scores.gap_extend - scores.gap_open];
                if (d - 1 >= wf_prev.diag_begin  && d - 1 < wf_prev.diag_begin + (int64_t) wf_prev.M.size()) {
                    if (wf.I[d - wf.diag_begin] == wf_prev.M[d - 1 - wf_prev.diag_begin] + 1) {
                        // an insert opened here
                        s -= scores.gap_extend + scores.gap_open;
                        op_len += 1;
                        cigar.emplace_back('I', op_len);
                        d -= 1;
                        a -= 1;
                        mat = MAT_M;
                        op_len = 0;
                        continue;
                    }
                }
            }
        }
        else {
            // we're in the deletion matrix
            if (s >= scores.gap_extend) {
                const auto& wf_prev = wfs[s - scores.gap_extend];
                if (d + 1 >= wf_prev.diag_begin && d + 1 < wf_prev.diag_begin + (int64_t) wf_prev.M.size()) {
                    if (wf.D[d - wf.diag_begin] == wf_prev.D[d + 1 - wf_prev.diag_begin] + 1) {
                        // a deletion extended here
                        s -= scores.gap_extend;
                        op_len += 1;
                        d += 1;
                        a -= 1;
                        continue;
                    }
                }
            }
            if (s >= scores.gap_extend + scores.gap_open) {
                const auto& wf_prev = wfs[s - scores.gap_extend - scores.gap_open];
                if (d + 1 >= wf_prev.diag_begin && d + 1 < wf_prev.diag_begin + (int64_t) wf_prev.M.size()){
                    if (wf.D[d - wf.diag_begin] == wf_prev.M[d + 1 - wf_prev.diag_begin] + 1) {
                        // a deletion opened here
                        s -= scores.gap_extend + scores.gap_open;
                        op_len += 1;
                        d += 1;
                        a -= 1;
                        cigar.emplace_back('D', op_len);
                        mat = MAT_M;
                        op_len = 0;
                        continue;
                    }
                }
            }
        }
    }

    // handle the final operation
    if (op_len != 0) {
        if (mat == MAT_M) {
            cigar.emplace_back('M', op_len);
        }
        else if (mat == MAT_D) {
            cigar.emplace_back('D', op_len);
        }
        else {
            cigar.emplace_back('I', op_len);
        }
    }

    std::reverse(cigar.begin(), cigar.end());
    if (cigar.back().len == 0) {
        cigar.pop_back();
    }
}

void wavefront_prune(Wavefront& wf, int32_t sub_opt_diff) {

    if (!wf.M.empty()) {
        int32_t max_anti_diag = *std::max_element(wf.M.begin(), wf.M.end());

        while (wf.M.front() < max_anti_diag - sub_opt_diff) {
            wf.M.pop_front();
            wf.I.pop_front();
            wf.D.pop_front();
            ++wf.diag_begin;
        }

        while (wf.M.back() < max_anti_diag - sub_opt_diff) {
            wf.M.pop_back();
            wf.I.pop_back();
            wf.D.pop_back();
        }
    }
}

void print_wfs(const std::vector<Wavefront>& wfs) {
    std::cerr << "# print wfs #" << std::endl;
    int min_begin = std::numeric_limits<int>::max();
    int max_end = std::numeric_limits<int>::min();
    for (int s = 0; s < wfs.size(); ++s) {
        min_begin = std::min<int>(wfs[s].diag_begin, min_begin);
        max_end = std::max<int>(wfs[s].diag_begin + wfs[s].M.size(), max_end);
    }
    for (int d = min_begin; d < max_end; ++d) {
        std::cerr << d << "\t";
    }
    std::cerr << std::endl;
    for (int s = 0; s < wfs.size(); ++s) {
        std::cerr << "score " << s << std::endl;
        for (int d = min_begin; d < max_end; ++d) {
            if (d >= wfs[s].diag_begin && d < wfs[s].diag_begin +  (int) wfs[s].M.size()) {
                int a = wfs[s].M[d - wfs[s].diag_begin];
                if (a > -100) {
                    std::cerr << a;
                }
                else {
                    std::cerr << ".";
                }
            }
            std::cerr << "\t";
        }
        std::cerr << std::endl;
        for (int d = min_begin; d < max_end; ++d) {
            if (d >= wfs[s].diag_begin && d < wfs[s].diag_begin +  (int) wfs[s].M.size()) {
                int a = wfs[s].I[d - wfs[s].diag_begin];
                if (a > -100) {
                    std::cerr << a;
                }
                else {
                    std::cerr << ".";
                }
            }
            std::cerr << "\t";
        }
        std::cerr << std::endl;
        for (int d = min_begin; d < max_end; ++d) {
            if (d >= wfs[s].diag_begin && d < wfs[s].diag_begin +  (int) wfs[s].M.size()) {
                int a = wfs[s].D[d - wfs[s].diag_begin];
                if (a > -100) {
                    std::cerr << a;
                }
                else {
                    std::cerr << ".";
                }
            }
            std::cerr << "\t";
        }
        std::cerr << std::endl;
    }
}

void wavefront_viz(const std::string& seq1, const std::string& seq2,
                   const std::vector<Wavefront>& wfs, path output_path) {
    size_t base_size = min(size_t(2000), max(seq1.size(), seq2.size()));

    size_t x_size;
    size_t y_size;

    if (seq1.size() > seq2.size()) {
        x_size = ceil(double(seq2.size())/seq1.size()) * base_size;
        y_size = base_size;
    }
    else{
        x_size = base_size;
        y_size = ceil(double(seq2.size())/seq1.size()) * base_size;
    }

    int32_t score_interval = wfs.size() / min(size_t(200), wfs.size());

    cerr << "x_size: " << x_size << '\n';
    cerr << "y_size: " << y_size << '\n';

    SvgPlot plot(output_path, x_size, y_size, 0, x_size, 0, y_size);
    string type = "rect";
    int16_t width = 1;

    Viridis viridis;

//    std::vector<std::vector<int>> matrix(seq1.size() + 1, std::vector<int32_t>(seq2.size() + 1, -1));

    vector <vector <bool> > is_filled(x_size, vector<bool>(y_size, false));

    for (int s = wfs.size() - 1; s >= 0; s -= score_interval) {
        auto color = viridis.get_svg_color(double(s)/double(wfs.size()));

        const auto& wf = wfs[s];
        for (int k = 0; k < wf.M.size(); ++k) {
            int d = wf.diag_begin + k;
            int a = wf.M[k];
            if (a < -2) {
                continue;
            }
            int i = (a + d) / 2;
            int j = (a - d) / 2;
            while (i >= 0 && j >= 0 && i < seq1.size() && j < seq2.size() && seq1[i] == seq2[j] && a != wf.I[k] && a != wf.D[k]) {
//                matrix[i + 1][j + 1] = s;

                size_t x = round(x_size*(double(i+1)/seq1.size()));
                size_t y = round(y_size*(double(j+1)/seq2.size()));

                if (x < x_size and y < y_size){
                    if (not is_filled[x][y]) {
                        plot.add_point(x, y, type, width, color);
                        is_filled[x][y] = true;
                    }
                }

                --i;
                --j;
                a -= 2;
            }
//            matrix[i + 1][j + 1] = s;

            size_t x = round(x_size*(double(i+1)/seq1.size()));
            size_t y = round(y_size*(double(j+1)/seq2.size()));

            if (x < x_size and y < y_size){
                if (not is_filled[x][y]) {
                    plot.add_point(x, y, type, width, color);
                    is_filled[x][y] = true;
                }
            }
        }
    }
}

int score_cigar(const std::string& seq1, const std::string& seq2,
                const std::vector<CIGAROp>& cigar, const WFScores& scores) {
    int s = 0;
    int i = 0;
    int j = 0;
    for (auto c : cigar) {
        switch (c.op) {
            case 'I':
                i += c.len;
                s += scores.gap_open + c.len * scores.gap_extend;
                break;
            case 'D':
                j += c.len;
                s += scores.gap_open + c.len * scores.gap_extend;
                break;
            case 'M':
            {
                for (int k = 0; k < c.len; ++k, ++i, ++j) {
                    if (seq1[i] != seq2[j]) {
                        s += scores.mismatch;
                    }
                }
            }            default:
                break;
        }
    }
    if (not (cigar.empty() || i == seq1.size())){
        throw runtime_error("ERROR: assertion failed: cigar.empty() || i == seq1.size()");
    }
    if (not (cigar.empty() || j == seq2.size())){
        throw runtime_error("ERROR: assertion failed: cigar.empty() || j == seq2.size()");
    }
    return s;
}

std::vector<CIGAROp> wavefront_align(const std::string& seq1, const std::string& seq2,
                                     const WFScores& scores, double max_score_a,
                                     int32_t max_score_b, int32_t prune_diff, path viz_path) {

    std::vector<CIGAROp> cigar;

    int32_t final_diag = seq1.size() - seq2.size();
    int32_t final_anti_diag = seq1.size() + seq2.size() - 2;

    // init wavefront to the upper left of both sequences' starts
    std::vector<Wavefront> wfs(1, Wavefront(0, 1));
    wfs.back().M[0] = -2;

    int32_t max_score = max_score_b;

    // do wavefront iterations until hit max score or alignment finishes
    wavefront_extend(seq1, seq2, wfs.back());
    while (!wavefront_reached(wfs.back(), final_diag, final_anti_diag)
           && wfs.size() <= max_score) {
        wavefront_next(seq1, seq2, scores, wfs);
        wavefront_extend(seq1, seq2, wfs.back());

        if (prune_diff >= 0) {
            // prune lagging diagonals
            wavefront_prune(wfs.back(), prune_diff);
        }
        if (max_score_a != 0.0 && !wfs.back().M.empty()) {
            // adjust the max score to the current
            int32_t max_anti_diag = std::max<int32_t>(0, *std::max_element(wfs.back().M.begin(), wfs.back().M.end()));
            max_score = max_score_a * max_anti_diag + max_score_b;
        }
    }

    if (not viz_path.empty()) {
        wavefront_viz(seq1, seq2, wfs, viz_path);
    }

    if (wfs.size() <= max_score) {
        // we didn't bail out for too high a score, do a traceback
        wavefront_traceback(seq1, seq2, scores, wfs, cigar);
    }

    return cigar;
}

