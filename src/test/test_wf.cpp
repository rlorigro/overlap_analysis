//
//  wf_test.cpp
//  
//
//  Created by Jordan Eizenga on 4/5/21.
//

#include <iostream>

#include "wavefront.hpp"

int main(void) {

    WFScores scores(1, 1, 1);

    // 0    2   3   4  6         7  8     11     12
    //          *   *            *  *            *
    // TACGGTCACCGCGACG-ACGACGGCAATTACTCCAAGTTGTCT
    // TACGG-CACGGCGATGAACGACGGCACTTTCTCCA--TTGTCC
    std::string seq1 = "TACGGTCACCGCGACGACGACGGCAATTACTCCAAGTTGTCT";
    std::string seq2 = "TACGGCACGGCGATGAACGACGGCACTTTCTCCATTGTCC";

    double a = 0.0;
    int b = 10;
    int prune = -1;
    path viz_path = "test.svg";
    std::vector<CIGAROp> cigar = wavefront_align(seq1, seq2, scores, a, b, prune, viz_path);

    std::cout << "sequence 1: " << std::endl;
    std::cout << seq1 << std::endl;
    std::cout << "sequence 2: " << std::endl;
    std::cout << seq2 << std::endl;
    std::cout << "CIGAR string:" << std::endl;
    for (auto cigar_op : cigar) {
        std::cout << cigar_op.len << cigar_op.op;
    }
    std::cout << std::endl;
}