#ifndef OVERLAP_ANALYSIS_SVGPLOT_HPP
#define OVERLAP_ANALYSIS_SVGPLOT_HPP

#include <experimental/filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::absolute;
using std::experimental::filesystem::exists;
using std::experimental::filesystem::path;
using std::runtime_error;
using std::ofstream;
using std::to_string;
using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::array;


class SvgPlot{
public:
    /// Attributes ///
    path file_path;
    path image_output_path;
    ofstream file;

    size_t width;
    size_t height;

    size_t x_min;
    size_t x_max;
    size_t y_min;
    size_t y_max;

    bool axes;
    double padding;

    /// Methods ///
    SvgPlot(
            path output_path,
            size_t width,
            size_t height,
            size_t x_min,
            size_t x_max,
            size_t y_min,
            size_t y_max,
            bool axes=false);

    ~SvgPlot();

    void check_file();

    template <class T, class T2, class T3> void add_line(
            const T x,
            const T y,
            const T x2,
            const T y2,
            T2 width,
            T3& color);

    template <class T, class T2, class T3> void add_lines(
            const vector <array <T, 2> >& coordinates,
            T2 width,
            T3& color);

    template <class T, class T2, class T3> void add_lines(
            const vector <T>& x,
            const vector <T>& y,
            T2 width,
            T3& color);

    template <class T, class T2, class T3> void add_disjoint_lines(
            const vector <array <T, 4> >& coordinates,
            T2 width,
            T3& color);

    template <class T, class T2, class T3> void add_points(
            const vector <array <T, 2> >& coordinates,
            string& type,
            T2 size,
            T3& color);

    template <class T, class T2, class T3> void add_point(
            T x,
            T y,
            string& type,
            T2 size,
            T3& color);
};


template <class T, class T2, class T3> void SvgPlot::add_line(
        const T x,
        const T y,
        const T x2,
        const T y2,
        T2 width,
        T3& color){

    check_file();

    file << '\t' << "<line x1='" << x << "' y1='" << y
         << "' x2='" << x2 << "' y2='" << y2
         << "' stroke='" << color << "' stroke-width='" << width << "' />" << '\n';
}


template <class T, class T2, class T3> void SvgPlot::add_lines(
        const vector <array <T, 2> >& coordinates,
        T2 width,
        T3& color){

    check_file();

    for (size_t i=0; i<coordinates.size() - 1; i++) {
        file << '\t' << "<line x1='" << coordinates[i][0] << "' y1='" << coordinates[i][1]
             << "' x2='" << coordinates[i+1][0] << "' y2='" << coordinates[i+1][1]
             << "' stroke='" << color << "' stroke-width='" << width << "' />" << '\n';
    }
}


template <class T, class T2, class T3> void SvgPlot::add_lines(
        const vector <T>& x,
        const vector <T>& y,
        T2 width,
        T3& color){

    check_file();

    if (x.size() != y.size()){
        throw runtime_error("ERROR: x and y vectors must be equal size");
    }

    for (size_t i=0; i<x.size() - 1; i++) {
        file << '\t' << "<line x1='" << x[i] << "' y1='" << y[i]
             << "' x2='" << x[i+1] << "' y2='" << y[i+1]
             << "' stroke='" << color << "' stroke-width='" << width << "' />" << '\n';
    }
}


template <class T, class T2, class T3> void SvgPlot::add_disjoint_lines(
        const vector <array <T, 4> >& coordinates,
        T2 width,
        T3& color){

    check_file();

    for (const auto& c: coordinates){
        file << '\t' << "<line x1='" << c[0] << "' y1='" << c[1] << "' x2='" << c[2] << "' y2='" << c[3]
             << "' stroke='" << color << "' stroke-width='" << width << "' />\n";
    }
}


template <class T, class T2, class T3> void SvgPlot::add_points(
        const vector <array <T, 2> >& coordinates,
        string& type,
        T2 size,
        T3& color){

    check_file();

    for (const auto& c: coordinates){
        double x = c[0] - double(size)/2;
        double y = c[1] - double(size)/2;

        file << '\t' << "<" << type << " x='" << x << "' y='" << y << "' width='" << size << "' height='" << size
             << "' stroke='" << color << "'/>" << '\n';
    }
}


template <class T, class T2, class T3> void SvgPlot::add_point(
        T x,
        T y,
        string& type,
        T2 size,
        T3& color){

    check_file();

    if (type == "circle") {
        file << '\t' << "<circle cx='" << x << "' cy='" << y << "' r='" << size
             << "' fill='" << color << "'/>" << '\n';
    }
    else if (type == "rect") {
        file << '\t' << "<rect x='" << x << "' y='" << y << "' width='" << size << "' height='" << size
             << "' fill='" << color << "'/>" << '\n';
    }
}


#endif //OVERLAP_ANALYSIS_PLOT_HPP
