#ifndef OVERLAP_ANALYSIS_FASTQITERATOR_HPP
#define OVERLAP_ANALYSIS_FASTQITERATOR_HPP

#include <experimental/filesystem>
#include <fstream>
#include <string>
#include <vector>

using std::experimental::filesystem::create_directories;
using std::experimental::filesystem::path;
using std::ifstream;
using std::to_string;
using std::string;
using std::vector;


double quality_char_to_error_probability(char q);


class FastqElement{
public:
    /// Attributes ///
    string name;
    string sequence;
    string quality_string;

    /// Methods ///
    FastqElement()=default;
    void get_error_probabilities(vector<float>& probabilities);
};


class FastqIterator {
public:
    /// Attributes ///
    const path file_path;

    /// Methods ///
    FastqIterator(path file_path);

    bool next_fastq_element(FastqElement& element);
    size_t get_line_number() const;

private:
    /// Attributes ///
    bool eof;
    ifstream file;
    size_t line_number;

};


#endif //OVERLAP_ANALYSIS_FASTQITERATOR_HPP
