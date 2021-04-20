#include "SvgPlot.hpp"


int main(){
    SvgPlot plot("test_lines.svg", 800, 600, 0, 6, 0, 6);

    vector<size_t> x = {0,5,6,4,3};
    vector<size_t> y = {5,1,4,6,3};

    string color = "blue";
    string color2 = "red";

    plot.add_lines(x, y, 0.25, color);

    SvgPlot plot_w_axes("test_lines_w_axes.svg", 800, 600, 0, 6, 0, 6, true);

    plot_w_axes.add_lines(x, y, 0.25, color);

    SvgPlot plot2("test_disjoint_lines.svg", 800, 600, 0, 6, 0, 6);

    vector <array <size_t,4> > coords = {
            {4,4,0,0},
            {4,4,0,2},
            {4,4,0,4},
            {4,4,2,0},
            {4,4,2,4},
            {4,4,4,0},
            {4,4,4,2},
            {4,4,4,4}
    };

    plot2.add_disjoint_lines(coords, 0.25, color);

    SvgPlot plot3("test_disjoint_curves.svg", 800, 600, 0, 6, 0, 6);

    plot3.add_disjoint_lines(coords, 0.025, color2);
    plot3.add_disjoint_curves(coords, 0.5, 0.05, color);

    return 0;
}
