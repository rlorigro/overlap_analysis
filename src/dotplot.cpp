#include "dotplot.hpp"
#include "SvgPlot.hpp"
#include "Plot.hpp"

#include "mummer/sparseSA.hpp"


void generate_dotplot(
        string& name,
        string& ref_sequence,
        string& query_sequence,
        uint32_t min_length,
        path& output_directory){

    path plot_path = output_directory / (name + "_dotplot.svg");
    SvgPlot plot(plot_path, 800, 800, 0, ref_sequence.size(), 0, ref_sequence.size());

    vector <array <int64_t,4> > coords;

    auto matcher = mummer::mummer::sparseSA::create_auto(ref_sequence.c_str(), ref_sequence.size(), 0, true);

    vector<mummer::mummer::match_t> mams;

    string type = "rect";
    string color = "#8115A3";

    matcher.findMAM_each(query_sequence, min_length, false, [&](const mummer::mummer::match_t& match){
        mams.emplace_back(match);

        for (size_t i=0; i<match.len; i+=50) {
            plot.add_point(
                    match.ref + i,
                    match.query + i,
                    type,
                    ref_sequence.size()/1000,
                    color);
        }
    });
}
