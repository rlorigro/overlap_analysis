#include "AlignmentChain.hpp"
#include <experimental/filesystem>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::path;



int main(){
//    path current_path = __FILE__;
//    path project_directory = current_path.parent_path().parent_path().parent_path();
//    path test_path = project_directory / "data" / "20_reads_VS_hg38_chr5_subset.sorted.paf";

    path test_path = "/home/ryan/data/nanopore/human/test/overlap_guppy_360_linear_subset_1/HG002_guppy360_id5001-1_radius40_subset_VS_HG002_mat_primary_sv-off.paf";

    AlignmentChains a;
    a.load_from_paf(test_path);

    a.split_all_chains();

    return 0;
}