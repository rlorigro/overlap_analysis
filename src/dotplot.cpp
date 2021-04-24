#include "dotplot.hpp"
#include "SvgPlot.hpp"
#include "Plot.hpp"

#include "mummer/sparseSA.hpp"



char complement(char c){
    c = char(toupper(c));

    if (c == 'A'){
        return 'T';
    }
    else if (c == 'C'){
        return 'G';
    }
    else if (c == 'G'){
        return 'C';
    }
    else if (c == 'T'){
        return 'A';
    }
    else {
        throw runtime_error("ERROR: non DNA character cannot be complemented: " + c);
    }
}


void reverse_complement(string& sequence, string& rc_sequence){
    rc_sequence.resize(0);

    for (std::string::reverse_iterator iter=sequence.rbegin(); iter!=sequence.rend(); ++iter){
        rc_sequence += complement(*iter);
    }
}


void generate_dotplot(
        const string& name,
        const string& ref_sequence,
        const string& query_sequence,
        uint32_t min_length,
        SvgPlot& plot,
        bool both_strands,
        int64_t x_offset,
        int64_t y_offset){

    vector <array <int64_t,4> > coords;

    auto matcher = mummer::mummer::sparseSA::create_auto(ref_sequence.c_str(), ref_sequence.size(), 0, true);

    vector<mummer::mummer::match_t> matches;

    string type = "rect";
    string color = "#8115A3";

    matcher.findMAM_each(query_sequence, min_length, false, [&](const mummer::mummer::match_t& match){
        matches.emplace_back(match);

        for (size_t i=0; i<match.len; i+=50) {
            plot.add_point(
                    match.ref + i + x_offset,
                    match.query + i + y_offset,
                    type,
                    ref_sequence.size()/1000,
                    color);
        }
    });

    if (both_strands){
        string rc_query;

    }
}
