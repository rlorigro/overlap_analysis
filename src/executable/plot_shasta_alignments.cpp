#include "SvgPlot.hpp"

#include "boost/program_options.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::bool_switch;
using boost::program_options::value;

#include <istream>
#include <iostream>
#include <functional>

using std::experimental::filesystem::directory_iterator;
using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::absolute;
using std::experimental::filesystem::rename;
using std::experimental::filesystem::exists;
using std::experimental::filesystem::path;
using std::ifstream;
using std::function;


class AlignmentCoordinate {
public:
    size_t coord_a;
    size_t coord_b;

    AlignmentCoordinate()=default;
};


void for_each_alignment_coordinate(path file_path, const function<void(const AlignmentCoordinate&)>& f){
    ifstream file(file_path);

    if (not (file.good() and file.is_open())){
        throw runtime_error("ERROR: file could not be read: " + file_path.string());
    }

    string token;
    AlignmentCoordinate alignment_coordinate;

    uint64_t n_delimiters = 0;
    uint64_t n_lines = 0;
    char c;

    // Skip first line
    file.ignore(9999, '\n');

    while (file.get(c)){
        if (c == ',') {
            if (n_delimiters == 3){
                alignment_coordinate.coord_a = stoi(token);
            }
            else if (n_delimiters == 4){
                alignment_coordinate.coord_b = stoi(token);
            }

            token.resize(0);
            n_delimiters++;
        }
        else if (c == '\n'){
            if (n_delimiters < 4 and n_lines > 0){
                throw runtime_error(
                        "ERROR: file provided does not contain sufficient tab delimiters to be shasta alignment: " +
                        to_string(n_lines));
            }

            f(alignment_coordinate);

            token.resize(0);
            n_delimiters = 0;
            n_lines++;
        }
        else {
            token += c;
        }
    }

}


void plot_shasta_alignment(path input_dir, path output_dir){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory already exists");
    }
    else{
        create_directories(output_dir);
    }

    string type = "rect";
    string color = "red";

    for(auto& file_iter: directory_iterator(input_dir)){
        const path& file_path = file_iter.path();
        path output_path = output_dir / file_path.filename();
        output_path.replace_extension(".svg");

        cerr << file_path << '\n';

        vector<AlignmentCoordinate> coordinates;

        for_each_alignment_coordinate(file_path, [&](const AlignmentCoordinate& a){
            coordinates.emplace_back(a);
        });

        size_t a_length = coordinates.back().coord_a;
        size_t b_length = coordinates.back().coord_b;

        size_t a_start = coordinates[0].coord_a;
        size_t b_start = coordinates[0].coord_b;

        cerr << a_length << ' ' << a_start << ' ' << b_length << ' ' << b_start << '\n';

        SvgPlot plot(output_path, 1200, 1200, 0, a_length-a_start, 0, b_length-b_start);

        size_t i = 0;
        for (auto& item: coordinates){
            plot.add_point(item.coord_a - a_start, item.coord_b - b_start, type, 20, color);

            if (i < 100){
                cerr << item.coord_a - a_start << ' ' << item.coord_b - b_start << ' ' << item.coord_a << ' ' << item.coord_b << '\n';
            }

            i++;
        }
    }
}


int main(int argc, char* argv[]){
    path alignment_directory;
    path output_directory;

    options_description options("Arguments:");

    options.add_options()
            ("input_dir",
             value<path>(&alignment_directory)
             ->required(),
             "Alignments TO BE EVALUATED: path of directory containing Shasta alignment details (one file per alignment)\n")

            ("output_dir",
             value<path>(&output_directory)
             ->required(),
             "Where to dump output SVG, CSV, PAF, etc")
             ;

    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cerr << options << "\n";
        return 0;
    }
    notify(vm);

    plot_shasta_alignment(
            alignment_directory,
            output_directory);

    return 0;
}
