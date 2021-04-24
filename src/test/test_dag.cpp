#include <DagAligner.hpp>

using std::experimental::filesystem::exists;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path output_directory = project_directory / "output/";

    if (not exists(output_directory)){
        create_directories(output_directory);
    }

    vector <pair <coord_t,size_t> > matches = {
            {{0+1,0+1}, 2},
            {{1+1,0+1}, 2},
            {{2+1,2+1}, 4},
            {{2+1,3+1}, 3},
            {{6+1,6+1}, 2}
    };

    Dag dag(matches, 9, 9, 2);

    cerr << dag << '\n';

    path output_path = output_directory/"dag.svg";
    cerr << "Writing SVG: " << output_path << '\n';

    dag.write_to_svg(output_path, false);

    dag.compute_alignment();

    path output_path2 = output_directory/"dag_aligned.svg";
    cerr << "Writing SVG: " << output_path2 << '\n';

    dag.write_to_svg(output_path2, false);

    return 0;
}
