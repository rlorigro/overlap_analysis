#ifndef OVERLAP_ANALYSIS_COLOR_HPP
#define OVERLAP_ANALYSIS_COLOR_HPP


#include <vector>
#include <array>
#include <string>

using std::vector;
using std::array;
using std::string;


class ColorMap {
    /// Methods ///
    virtual array<double,3> get_rgb(double x)=0;
};


class Viridis: public ColorMap {
public:
    /// Attributes ///
    static const vector <array <double,3> > rgbs;

    /// Methods ///
    Viridis()=default;
    array<double,3> get_rgb(double x);
    string get_svg_color(double x);
};





#endif //OVERLAP_ANALYSIS_COLOR_HPP
