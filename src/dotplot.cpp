#include "dotplot.hpp"
#include "Plot.hpp"

#include "mummer/sparseSA.hpp"


void generate_dotplot(
        string& name,
        string& ref_sequence,
        string& query_sequence,
        uint32_t min_length,
        path& output_directory){

    path plot_path = output_directory / (name + "_dotplot.svg");
    Plot plot(plot_path, 800, 800, 18, 2);

    vector <array <int64_t,2> > coords;

    auto matcher = mummer::mummer::sparseSA::create_auto(ref_sequence.c_str(), ref_sequence.size(), 0, true);

    vector<mummer::mummer::match_t> mams;

    matcher.findMAM_each(query_sequence, min_length, false, [&](const mummer::mummer::match_t& match){
        mams.emplace_back(match);

        for (int64_t i=0; i<match.len; i+=50){
            array <int64_t,2> coord = {match.ref + i, match.query + i};

            coords.emplace_back(coord);
        }
    });

    plot.add_points(coords, "", 0.1, "");
    plot.generate();
}
