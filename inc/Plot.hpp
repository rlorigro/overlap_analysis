#ifndef OVERLAP_ANALYSIS_PLOT_HPP
#define OVERLAP_ANALYSIS_PLOT_HPP

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


class Plot{
public:
    path file_path;
    path directory = "/dev/shm/";
    const string filename_suffix = ".gp";
    ofstream file;

    Plot(path output_path, size_t width, size_t height);
    void check_file();

    void set_xrange(size_t x_min, size_t x_max);
    void set_yrange(size_t y_min, size_t y_max);

    void generate();
    template <class T> void add_lines(const vector<T>& x, const vector<T>& y, string line_color="blue", string title="");
    template <class T> void add_disjoint_lines(const vector <array <T, 4> >& coordinates, string line_color="blue", string title="");
};


template <class T> void Plot::add_lines(const vector<T>& x, const vector<T>& y, string line_color, string title) {
    check_file();

    if (x.size() != y.size()){
        throw runtime_error("ERROR: can't plot unequal size x and y vectors: " + to_string(x.size()) + " " + to_string(y.size()));
    }

    file << "plot '-' with lines ";

    if (not title.empty()) {
        file << "title " << title << '\n';
    }
    else{
        file << "notitle" << '\n';
    }

    for (size_t i=0; i<x.size(); i++){
        file << '\t' << x[i] << "," << y[i] << '\n';
    }

    // Denote end of data
    file << "\te\n";
}


template <class T> void Plot::add_disjoint_lines(const vector <array <T, 4> >& coordinates, string line_color, string title) {
    check_file();

    file << "plot '-' u 1:2:3:4 w vectors nohead ";

    if (not title.empty()) {
        file << "title " << title << '\n';
    }
    else{
        file << "notitle" << '\n';
    }

    for (const auto& c: coordinates){
        file << '\t' << c[0] << ' ' << c[1] << ' ' << c[2] - c[0] << ' ' << c[3] - c[1] << '\n';
    }

    // Denote end of data
    file << "\te\n";
}



#endif //OVERLAP_ANALYSIS_PLOT_HPP
