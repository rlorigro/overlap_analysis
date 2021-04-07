#include "Plot.hpp"

#include "boost/uuid/random_generator.hpp"
#include "boost/uuid/uuid_io.hpp"
#include "boost/lexical_cast.hpp"

using boost::uuids::random_generator;
using std::to_string;


void Plot::check_file(){
    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: cannot plot with closed gnuplot file");
    }
}


Plot::Plot(path output_path, size_t width, size_t height){
    // Maybe there is a collision with another thread, so allow several retries on generating uuid
    int8_t retries = 5;
    while (--retries >= 0){
        if (file.is_open() and file.good()){
            break;
        }
        else{
            file_path = directory / (to_string(random_generator()()) + filename_suffix);
            file.open(file_path);
        }
    }

    file << "set terminal svg size " << width << "," << height << " font 'Noto Serif'\n"
         << "set output " << absolute(output_path) << "\n"
         << "set xtics out nomirror\n"
         << "set ytics out nomirror\n";
}


void Plot::set_xrange(size_t x_min, size_t x_max){
    check_file();

    file << "set xrange [" << x_min << ',' << x_max << "]\n";
}


void Plot::set_yrange(size_t y_min, size_t y_max){
    check_file();

    file << "set yrange [" << y_min << ',' << y_max << "]\n";
}


void Plot::add_lines(const vector<size_t>& x, const vector<size_t>& y, string line_color, string title) {
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


void Plot::generate(){
    check_file();

    file.close();

    const string command = "gnuplot " + file_path.string();

    ::system(("echo " + command).c_str());
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("ERROR " +
                            to_string(errorCode) + " " + strerror(errorCode) +
                            "\nrunning command: " + command);
    }
}
