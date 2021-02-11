#include "FastqIterator.hpp"
#include <iostream>
#include <stdexcept>

using std::runtime_error;
using std::cout;
using std::cerr;

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

using namespace boost::accumulators;

typedef accumulator_set<float, features<tag::count, tag::mean, tag::variance>> stats_accumulator;


void split_palindrome_by_quality(FastqElement& element, vector<FastqElement>& result){
    stats_accumulator left_stats;
    stats_accumulator right_stats;

    auto length = element.quality_string.size();

    if (length < 6){
        return;
    }

    size_t midpoint = length/2;

    for (size_t i=0; i<midpoint; i++){
        auto q = element.quality_string[i];
        auto p = quality_char_to_error_probability(q);

        left_stats(p);
    }

    for (size_t i=midpoint; i<length; i++){
        auto q = element.quality_string[i];
        auto p = quality_char_to_error_probability(q);

        right_stats(p);
    }

    float left_mean = mean(left_stats);
    float left_variance = variance(left_stats);

    float right_mean = mean(right_stats);
    float right_variance = variance(right_stats);

//    cout << left_mean << '\t' << left_variance << '\n';
//    cout << right_mean << '\t' << right_variance << '\n';

    if (right_mean - left_mean > 0.09 and right_mean >= 0.15){
        if (right_variance > left_variance and right_variance > 0.025){
            cout << '@' << element.name << '\n';
            cout << element.sequence << '\n';
            cout << '+' << '\n';
            cout << element.quality_string << '\n';
        }
    }
}


int main(int argc, char **argv){
    if (argc == 1 or argc > 2){
        cerr << "Usage: split_palindromes_by_quality /path/to/file.fastq\n";
        return 1;
    }

    path fastq_path = argv[1];
    FastqIterator iterator(fastq_path);
    FastqElement e;
    vector<FastqElement> result;

//    size_t i = 0;
    while (iterator.next_fastq_element(e)){
//        cout << ++i << " " << e.name << '\n';

        split_palindrome_by_quality(e, result);
    }

    return 0;
}
