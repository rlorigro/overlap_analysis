#include "boost/program_options.hpp"
#include "boost/icl/interval_map.hpp"
#include "boost/icl/interval.hpp"
#include "boost/bimap.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using boost::icl::total_enricher;
using boost::icl::interval_map;
using boost::icl::interval;
using boost::bimap;

#include <experimental/filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

using std::experimental::filesystem::path;
using std::ifstream;
using std::string;
using std::cerr;
using std::cout;
using std::map;


typedef map <string, interval_map<uint64_t,bool,total_enricher> > regional_interval_map;
typedef bimap<uint32_t,uint32_t> uint32_bimap;
typedef uint32_bimap::value_type bimap_pair;


void get_read_id_map(path& read_csv_path, uint32_bimap& id_vs_name){
    ///
    /// Parse a line of the ReadSummary.csv with the following format:
    ///     Id,Name, ... etc
    ///     0,512, ... etc
    ///
    ifstream read_info_file(read_csv_path);

    string token;
    uint32_t id;
    uint32_t name;
    uint64_t n_commas = 0;
    uint64_t n_lines = 0;
    char c;

    while (read_info_file.get(c)) {
        if (c == ',') {
            if (n_lines == 0){
                continue;
            }

            if (n_commas == 0) {
                id = stoi(token);
            }
            else if (n_commas == 1) {
                name = stoi(token);
                id_vs_name.insert(bimap_pair(id, name));
            }

            token.resize(0);
            n_commas++;
        }
        else if (c == '\n'){
            token.resize(0);
            n_commas = 0;
            n_lines++;
        }
        else {
            token += c;
        }
    }
}


void label_alignment_info(path info_csv_path, path read_csv_path){
    uint32_bimap id_vs_name;
    get_read_id_map(read_csv_path, id_vs_name);

    cerr << 37 << "->" << id_vs_name.left.at(37) << '\n';
}


int main(int argc, char* argv[]){
    path info_csv_path;
    path read_csv_path;

    options_description options("Arguments:");

    options.add_options()
            ("info_csv_path",
             value<path>(&info_csv_path),
             "File path of FASTA file containing query sequences to be aligned")

            ("read_csv_path",
             value<path>(&read_csv_path),
             "File path of FASTA file containing reference sequences to be aligned");

    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cerr << options << "\n";
        return 0;
    }

    label_alignment_info(info_csv_path, read_csv_path);

    return 0;
}
