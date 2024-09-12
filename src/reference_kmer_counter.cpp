#include "kmc_file.h"
#include "kmer_api.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <zlib.h>
#include <getopt.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
KSEQ_INIT(gzFile, gzread)

void usage() {
    std::cout << " Usage: ./reference_kmer_counter <args>\n";
    std::cout << " Mandatory args:\n";
    std::cout << "  -i, --kmc_db <str> \t \t KMC database path (without pre/suf).\n";
    std::cout << "  -r, --reference <str> \t \t The fasta file of the reference to scan (accepts gzipped file).\n";
    std::cout << " Optional args:\n";
    std::cout << "  -o, --output_prefix <str> \t Output prefix: output will be PREFIXkmer_count.tsv and PREFIXkmer_histo.tsv [DEAFULT: ""]\n";
    std::cout << "  -b, --bin_size <int> \t \t Bin size [DEAFULT: 10000]\n";
    std::cout << "  -s, --bin_shift <int> \t Bin sliding window shift [DEAFULT: 10000]\n";
    std::cout << "  -m, --min_kmer <int> \t \t Min kmer count [DEAFULT: 2]\n";
    std::cout << "  -M, --max_kmer <int> \t \t Max kmer count [DEAFULT: 200]\n";
    std::cout << "  -h, --help \t \t \t Prints the usage.\n";
}

struct InputFlags {
    std::string kmc_database_path;
    std::string reference_path;
    std::string statistic_output_path;
    std::string histogram_output_path;
    size_t bin_size;
    size_t bin_shift_size;
    uint64_t min_kmer_counter;
    uint64_t max_kmer_counter;
};

static struct option long_options[] = {
    {"kmc_db", required_argument, NULL, 'i'},
    {"reference", required_argument, NULL, 'r'},
    {"output_prefix", required_argument, NULL, 'o'},
    {"bin_size", required_argument, NULL, 'b'},
    {"bin_shift", required_argument, NULL, 's'},
    {"min_kmer", required_argument, NULL, 'm'},
    {"max_kmer", required_argument, NULL, 'M'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}};
    
void decodeFlags(int argc, char *argv[], InputFlags &flags) {
    int args = 0;
    int opt;
    int given_k = 0;
    
    std::string output_prefix = "";
    
    flags.bin_size = 10000;
    flags.bin_shift_size = 10000;
    
    flags.min_kmer_counter = 2;
    flags.max_kmer_counter = 200;
    
    bool is_kmc_path = false;
    bool is_ref_path = false;
    
    /* initialisation */
    while ((opt = getopt_long(argc, argv, "i:r:o:b:s:m:M:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'i':
                flags.kmc_database_path = optarg;
                is_kmc_path = true;
                args++;
                break;
            case 'r':
                flags.reference_path = optarg;
                is_ref_path = true;
                args++;
                break;
            case 'o':
                output_prefix = std::string(optarg);
                args++;
                break;
            case 'b':
                flags.bin_size = (size_t)atoi(optarg);
                args++;
                break;
            case 's':
                flags.bin_shift_size = (size_t)atoi(optarg);
                args++;
                break;
            case 'm':
                flags.min_kmer_counter = (size_t)atoi(optarg);
                args++;
                break;
            case 'M':
                flags.max_kmer_counter = (size_t)atoi(optarg);
                args++;
                break;
            case 'h':
                usage();
                exit(0);
            default:
                std::cerr << "Unknown option " << opt << std::endl;
                usage();
                exit(1);
        }
    }
    
    
    
    if(!is_kmc_path) {
        std::cerr << "Argument --kmc_db or -i needed." << std::endl;
        exit(1);
    }
    if(!is_ref_path) {
        std::cerr << "Argument --reference or -r needed." << std::endl;
        exit(1);
    }
    
    
    flags.statistic_output_path = output_prefix + "kmer_count.tsv";
    flags.histogram_output_path = output_prefix + "kmer_histo.tsv";
    
    
    std::string kmc_database_path;
    std::string reference_path;
    std::string statistic_output_path;
    std::string histogram_output_path;
    size_t bin_size;
    size_t bin_shift_size;
    uint64_t min_kmer_counter;
    uint64_t max_kmer_counter;
    
    std::cout << "KMC database path: " << flags.kmc_database_path << std::endl;
    std::cout << "Reference path: " << flags.reference_path << std::endl;
    
    std::cout << "Statistics output path: " << flags.statistic_output_path << std::endl;
    std::cout << "Histogram output path: " << flags.histogram_output_path << std::endl;
    
    std::cout << "Bin size: " << flags.bin_size << std::endl;
    std::cout << "Bin shift size: " << flags.bin_shift_size << std::endl;
    
    std::cout << "Min kmer counter: " << flags.min_kmer_counter << std::endl;
    std::cout << "Max kmer counter: " << flags.max_kmer_counter << std::endl;
}

int main(int argc, char **argv) {
    InputFlags flags;
    decodeFlags(argc, argv, flags);
    
    std::cout << "Reading database from: " << flags.kmc_database_path << std::endl;
    CKMCFile kmc_database;
    if(!kmc_database.OpenForRA(flags.kmc_database_path)) {
        std::cerr << "Failed to open KMC database." << std::endl;
        return 1;
    }
    size_t bin_size = flags.bin_size;
    size_t bin_shift_size = flags.bin_shift_size;
    
    uint64_t min_counter = flags.min_kmer_counter;
    uint64_t max_counter = flags.max_kmer_counter;
    
    CKMCFileInfo kmc_info;
    kmc_database.Info(kmc_info);
    uint32_t k = kmc_info.kmer_length;
    
    std::cout << "KMC database read. k=" << kmc_info.kmer_length << std::endl;
    
    CKmerAPI * kmer_object = new CKmerAPI(kmc_info.kmer_length);
    
    gzFile fp = gzopen(flags.reference_path.c_str(), "r");
    kseq_t *seq;
    seq = kseq_init(fp);
    int l;
    
    std::ofstream output;
    output.open(flags.statistic_output_path);
    
    std::ofstream histogram_output;
    histogram_output.open(flags.histogram_output_path);
    
    while((l = kseq_read(seq)) >= 0) {
        std::cerr << "Processing: " << seq->name.s << std::endl;
        std::string current_sequence(seq->seq.s);
        
        
        for(size_t i = 0; i < seq->seq.l; i += bin_shift_size) {
            uint64_t current_bin_total = 0;
            uint64_t count_below_min = 0;
            uint64_t count_above_max = 0;
            uint64_t count_n_kmer = 0;
            
            size_t max_idx = i + bin_size;
            
            int32_t last_n_pos = -1;
            
            std::vector<uint64_t> current_frequencies;
            std::vector<uint64_t> freq_histogram(max_counter - min_counter + 1);
            
            if(max_idx >= seq->seq.l) max_idx = seq->seq.l;
            for(size_t j = i; j < max_idx; j++) {
                if(current_sequence[j] == 'N' || current_sequence[j] == 'n') last_n_pos = j;
                if(j - i >= k && (last_n_pos == -1 || j - last_n_pos >= k)) {
                    kmer_object->from_string_impl(current_sequence.begin() + j, k);
                    uint64_t counter;
                    if(kmc_database.CheckKmer(*kmer_object, counter)) {
                        if(counter >= min_counter && counter <= max_counter) {
                            current_bin_total += counter;
                            current_frequencies.push_back(counter);
                        } else if(counter < min_counter) count_below_min++;
                        else count_above_max++;
                    } else {
                        kmer_object->reverse();
                        if(kmc_database.CheckKmer(*kmer_object, counter)) {
                            if(counter >= min_counter && counter <= max_counter) {
                                current_bin_total += counter;
                                current_frequencies.push_back(counter);
                            } else if(counter < min_counter) count_below_min++;
                            else count_above_max++;
                        }  else {
                            if(min_counter == 0) current_frequencies.push_back(0);
                        }
                    }
                }
                if(j - last_n_pos < k) count_n_kmer += 1;
            }
            
            if(current_frequencies.size() > 0) {
                std::sort(current_frequencies.begin(), current_frequencies.end());
                
                uint64_t current_number = current_frequencies[0];
                size_t current_count = 1;
                uint64_t current_mode = current_number;
                size_t current_highest = current_count;
                for(size_t j = 1; j < current_frequencies.size(); j++) {
                    if(current_frequencies[j] == current_number) current_count++;
                    else {
                        if(current_count > current_highest) {
                            current_highest = current_count;
                            current_mode = current_number;
                        }
                        current_number = current_frequencies[j];
                        current_count = 1;
                    }
                }
                if(current_count > current_highest) {
                    current_highest = current_count;
                    current_mode = current_number;
                }
                
                for(size_t j = 0; j < current_frequencies.size(); j++) freq_histogram[current_frequencies[j] - min_counter]++;
                
                uint64_t mode = current_mode;
                uint64_t median;
                if(current_frequencies.size() & 1) median = current_frequencies[current_frequencies.size() / 2];
                else if(current_frequencies.size() > 0) median = (current_frequencies[current_frequencies.size() / 2] + current_frequencies[current_frequencies.size() / 2 - 1]) / 2;
                
                output << seq->name.s << "\t" << i << "\t" << max_idx << "\t" << current_bin_total << "\t" << current_frequencies.size() << "\t" << double_t(current_bin_total)/double_t(current_frequencies.size()) << "\t" << mode << "\t" << median << "\t" << count_n_kmer << "\t" << count_below_min << "\t" << count_above_max << "\n";
                histogram_output << seq->name.s << "\t" << i << "\t" << max_idx;
                for(uint64_t j = 0; j < freq_histogram.size(); j++) histogram_output << "\t" << freq_histogram[j];
                histogram_output << "\n";
            } else {
                uint64_t num_bases = max_idx - i;
                output << seq->name.s << "\t" << i << "\t" << max_idx << "\t" << 0 << "\t" << 0 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << count_n_kmer << "\t" << count_below_min << "\t" << count_above_max << "\n";
                
                histogram_output << seq->name.s << "\t" << i << "\t" << max_idx;
                for(uint64_t j = 0; j < freq_histogram.size(); j++) histogram_output << "\t" << freq_histogram[j];
                histogram_output << "\n";
            }
        }
    }
    
    output.close();
    histogram_output.close();
    return 0;
}
