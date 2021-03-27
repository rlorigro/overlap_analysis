#include "AlignmentChain.hpp"
#include <experimental/filesystem>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::path;



int main(){
    path current_path = __FILE__;
    path project_directory = current_path.parent_path().parent_path().parent_path();
    path test_path = project_directory / "data" / "test_chain_long.paf";

    AlignmentChains a;
    a.load_from_paf(test_path);

    a.split_all_chains();

    return 0;
}