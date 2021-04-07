#ifndef OVERLAP_ANALYSIS_PLOT_HPP
#define OVERLAP_ANALYSIS_PLOT_HPP

#include <experimental/filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::absolute;
using std::experimental::filesystem::exists;
using std::experimental::filesystem::path;
using std::runtime_error;
using std::ofstream;
using std::string;
using std::vector;
using std::cerr;
using std::cout;


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
    void add_lines(const vector<size_t>& x, const vector<size_t>& y, string line_color="blue", string title="");
};


#endif //OVERLAP_ANALYSIS_PLOT_HPP
