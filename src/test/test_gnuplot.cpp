#include "Plot.hpp"


int main(){
    Plot plot("test.svg", 800, 600);

    vector<size_t> x = {1,5,2,4,3};
    vector<size_t> y = {5,1,4,2,3};

    plot.add_lines(x, y);
    plot.generate();

    return 0;
}