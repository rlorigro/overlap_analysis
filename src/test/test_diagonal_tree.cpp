#include "DiagonalTree.hpp"

#include <chrono>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <utility>

using std::vector;
using std::cerr;
using std::to_string;
using std::make_pair;
using std::pair;
using std::runtime_error;


void test(){
    vector <pair <pair <size_t,size_t>, size_t> > matches = {
            {{0,0}, 1},
            {{0,1}, 1},
            {{2,2}, 1},
            {{2,3}, 1},
            {{3,4}, 1},
            {{4,1}, 1},
    };

    DiagonalTree tree(5,5);

    for (const auto& m: matches){
        auto x = m.first.first;
        auto y = m.first.second;

        tree.insert(x, y, m.second);

        // Check that square and diagonal coordinates are interconvertible
        auto x_diag = tree.get_x_diagonal(x,y);
        auto y_diag = tree.get_y_diagonal(x,y);
        auto x1 = tree.get_x(x_diag,y_diag);
        auto y1 = tree.get_y(x_diag,y_diag);

        cerr << x << ' ' << y << ' ' << x_diag << ' ' << y_diag << ' ' << x1 << ' ' << y1 << '\n';
        if (not (x == x1 and y == y1)){
            throw runtime_error("ERROR: diagonal conversion failed for " + to_string(x) + ',' + to_string(y));
        }
    }

    size_t x_diag = 0;
    for (const auto& t: tree.diagonals){
        for(const auto& item: t){
            cerr << x_diag << ' ' << item.first << ' ' << item.second << '\n';
        }

        x_diag++;
    }

    for (const auto& m: matches){
        vector <pair <pair <size_t,size_t>, size_t> > results;

        cerr << "finding: " << m.first.first+m.second << ' ' << m.first.second+m.second << '\n';
        tree.find(m.first.first+m.second, m.first.second+m.second, 3, results);
        cerr << '\n';
    }
}


void time(){
    DiagonalTree tree(800,80);


    auto t1 = std::chrono::high_resolution_clock::now();

    for (size_t i=0; i<800; i++){
        for (size_t j=0; j<80; j++) {
            tree.insert(i,j,i*j);
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    cerr << "Construction completed in: " << elapsed.count() << " milliseconds" << '\n';

    vector<pair <pair <size_t, size_t>, size_t> > result;
    for (size_t i=0; i<800; i++){
        for (size_t j=0; j<80; j++) {
            tree.find(i,j,10,result,true);
        }
    }

    auto t3 = std::chrono::high_resolution_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
    cerr << "Queries completed in: " << elapsed2.count() << " milliseconds" << '\n';
}



int main()
{
    test();
//    time();
    return 0;
}
