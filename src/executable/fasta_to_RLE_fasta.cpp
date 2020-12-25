#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "SequenceElement.hpp"
#include "RunlengthSequenceElement.hpp"

#include <vector>
#include <thread>
#include <string>
#include <mutex>
#include <exception>
#include <atomic>
#include <iostream>
#include <experimental/filesystem>

#include "boost/program_options.hpp"

using std::vector;
using std::string;
using std::thread;
using std::cout;
using std::cerr;
using std::flush;
using std::mutex;
using std::lock_guard;
using std::thread;
using std::ref;
using std::move;
using std::exception;
using std::atomic;
using std::atomic_fetch_add;
using std::max;
using std::min;
using std::experimental::filesystem::path;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


template<class T> void runlength_encode(RunlengthSequenceElement& runlength_sequence, T& sequence){
    runlength_sequence = {};

    // First just copy the name
    runlength_sequence.name = sequence.name;

    char current_character = 0;

    // Iterate each run and append/update base/length for their corresponding vectors
    for (auto& character: sequence.sequence){
        if (tolower(character) != tolower(current_character)){
            runlength_sequence.sequence += character;
            runlength_sequence.lengths.push_back(1);
        }
        else{
            runlength_sequence.lengths.back()++;
        }

        current_character = character;
    }
}


template<class T> void runlength_decode(RunlengthSequenceElement& runlength_sequence, T& sequence){
    sequence = {};

    // First just copy the name
    runlength_sequence.name = sequence.name;

    // Iterate each run and append/update base/length for their corresponding vectors
    for (size_t i=0; i<runlength_sequence.sequence.size(); i++){
        sequence.sequence += string(runlength_sequence.lengths[i], runlength_sequence.sequence[i]);
    }
}


template<class T> void runlength_decode(RunlengthSequenceElement& runlength_sequence, T& sequence, size_t start, size_t stop){
    sequence = {};

    // First just copy the name
    runlength_sequence.name = sequence.name;

    // Iterate each run and append/update base/length for their corresponding vectors
    for (size_t i=start; i<stop; i++){
        sequence.sequence += string(runlength_sequence.lengths[i], runlength_sequence.sequence[i]);
    }
}


void runlength_encode_fasta_sequence_to_file(path& fasta_path,
                                             vector <pair <string, FastaIndex> >& read_index_vector,
                                             unordered_map<string, RunlengthSequenceElement>& runlength_sequences,
                                             mutex& map_mutex,
                                             mutex& file_write_mutex,
                                             FastaWriter& fasta_writer,
                                             bool store_in_memory,
                                             atomic<uint64_t>& job_index){

    while (job_index < read_index_vector.size()) {
        // Fetch add
        uint64_t thread_job_index = job_index.fetch_add(1);

        // Initialize containers
        SequenceElement sequence;
        RunlengthSequenceElement runlength_sequence;

        // Fetch Fasta sequence
        FastaReader fasta_reader = FastaReader(fasta_path);

        // Get read name and index
        string read_name = read_index_vector[thread_job_index].first;
        uint64_t read_index = read_index_vector[thread_job_index].second.byte_index;

        // First of each element is read_name, second is its index
        fasta_reader.get_sequence(sequence, read_name, read_index);

        // Convert to Run-length Encoded sequence element
        runlength_encode(runlength_sequence, sequence);

        // Write RLE sequence to file (no lengths written)
        file_write_mutex.lock();
        fasta_writer.write(runlength_sequence);
        file_write_mutex.unlock();

        if (store_in_memory) {
            // Append the sequence to a map of names:sequence
            map_mutex.lock();
            runlength_sequences[sequence.name] = move(runlength_sequence);
            map_mutex.unlock();
        }
        // Print status update to stdout
        cerr << "\33[2K\rParsed: " << sequence.name << flush;
    }
}


path runlength_encode_fasta_file(path input_file_path,
                                 unordered_map <string, RunlengthSequenceElement>& runlength_sequences,
                                 path output_dir,
                                 bool store_in_memory,
                                 uint16_t max_threads) {

    // Generate parent directories if necessary
    create_directories(output_dir);

    string output_filename;
    output_filename = string(input_file_path.filename());
    output_filename = output_filename.substr(0, output_filename.find_last_of(".")) + "_RLE.fasta";

    path output_file_path = output_dir / output_filename;

    cerr << "READING FILE: " << input_file_path.string() << "\n";
    cerr << "WRITING FILE: " << output_file_path.string() << "\n";

    // This reader is used to fetch an index
    FastaReader fasta_reader(input_file_path);

    // This writer is mutexed across threads
    FastaWriter fasta_writer(output_file_path);

    // Get index
    unordered_map <string, FastaIndex> read_indexes;
    fasta_reader.get_indexes_mapped_by_name(read_indexes);

    // Convert the map object into an indexable object
    vector <pair <string, FastaIndex> > read_index_vector;
    get_vector_from_index_map(read_index_vector, read_indexes);

    mutex map_mutex;
    mutex file_write_mutex;

    atomic<uint64_t> job_index = 0;
    vector<thread> threads;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        // Get data to send to threads (must not be sent by reference, or will lead to concurrency issues)
        try {
            // Call thread safe function to RL encode and write to file
            threads.emplace_back(thread(runlength_encode_fasta_sequence_to_file,
                                        ref(input_file_path),
                                        ref(read_index_vector),
                                        ref(runlength_sequences),
                                        ref(map_mutex),
                                        ref(file_write_mutex),
                                        ref(fasta_writer),
                                        store_in_memory,
                                        ref(job_index)));
        } catch (const exception &e) {
            cerr << e.what() << "\n";
            exit(1);
        }
    }

    // Wait for threads to finish
    for (auto& t: threads){
        t.join();
    }

    cerr << "\n" << flush;

    return output_file_path;
}


void fasta_to_RLE_fasta(path input_file_path, path output_dir, uint max_threads) {
    // Generate parent directories if necessary
    create_directories(output_dir);

    // Runlength encode the READ SEQUENCES, rewrite to another FASTA, and DON'T store in memory
    unordered_map<string, RunlengthSequenceElement> _;
    path reads_fasta_path_rle;
    bool store_in_memory = false;
    reads_fasta_path_rle = runlength_encode_fasta_file(
            input_file_path,
            _,
            output_dir,
            store_in_memory,
            max_threads);

    cout << '\n';
}


int main(int argc, char* argv[]){
    path input_file_path;
    path output_dir;
    uint16_t max_threads;

    options_description options("Required options");

    options.add_options()
            ("fasta",
             value<path>(&input_file_path),
             "File path of FASTA file containing sequences to be Run-length encoded")

            ("output_dir",
             value<path>(&output_dir)->
                     default_value("output/"),
             "Destination directory. File will be named based on input file name")

            ("max_threads",
             value<uint16_t>(&max_threads)->
                     default_value(1),
             "Maximum number of threads to launch");

    // Store options in a map and apply values to each corresponding variable
    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cout << options << "\n";
        return 0;
    }

    cout << "READING FILE: " << string(input_file_path) << "\n";

    fasta_to_RLE_fasta(input_file_path, output_dir, max_threads);

    return 0;
}

