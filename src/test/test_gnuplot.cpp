#include "Plot.hpp"


int main(){
    Plot plot("test_lines.svg", 800, 600);

    vector<size_t> x = {1,5,2,4,3};
    vector<size_t> y = {5,1,4,2,3};

    plot.add_lines(x, y);
    plot.generate();

    Plot plot2("test_disjoint_lines.svg", 800, 600);

    vector <array <size_t,4> > coords = {
            {0,0,1,1},
            {0,1,1,2},
            {0,2,1,3}
    };

    plot2.add_disjoint_lines(coords);
    plot2.generate();

    return 0;
}
