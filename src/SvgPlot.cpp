#include "SvgPlot.hpp"


SvgPlot::SvgPlot(path output_path, size_t width, size_t height, size_t x_min, size_t x_max, size_t y_min, size_t y_max):
        image_output_path(output_path),
        width(width),
        height(height),
        x_min(x_min),
        x_max(x_max),
        y_min(y_min),
        y_max(y_max),
        file(output_path)
{
    file << "<svg width='" << width << "' height='" << width << "' xmlns='http://www.w3.org/2000/svg' "
         << "viewBox='" << x_min << ' ' << y_min << ' ' << x_max << ' ' << y_max << "'>" << '\n';
}


SvgPlot::~SvgPlot(){
    file << "</svg>" << '\n';
    file.close();
}


void SvgPlot::check_file(){
    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: cannot plot with closed svg file");
    }
}
