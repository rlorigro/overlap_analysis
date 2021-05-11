#include "SweepLineSolver.hpp"

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
    size_t x_size = 10;
    size_t y_size = 10;

    vector <DiagonalMatch> matches = {
            {0,0,2,y_size},
            {2,1,2,y_size},
            {3,3,2,y_size},
            {3,0,2,y_size},
            {5,8,1,y_size},
    };

    for (const auto& item: matches){
        cerr << item.get_x(y_size) << " " << item.get_y(y_size) << '\n';
    }

    SweepLineSolver solver(matches);


}


void time(){
}



int main()
{
    test();
    return 0;
}
