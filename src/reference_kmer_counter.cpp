#include "kmc_file.h"
#include "kmer_api.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <zlib.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
KSEQ_INIT(gzFile, gzread)

int main(int argc, char **argv) {   
    if(argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <KMC database path> <Reference path> <Output file> <Bin size> <min kmer count> <max kmer count>" << std::endl;
        return 1;
    }
    
    std::cout << "Reading database from: " << argv[1] << std::endl;
    CKMCFile kmc_database;
    if(!kmc_database.OpenForRA(argv[1])) {
        std::cerr << "Failed to open KMC database." << std::endl;
        return 1;
    }
    
    std::cout << "Reading sequences from: " << argv[2] << std::endl;
    std::cout << "Outputting to: " << argv[3] << std::endl;
    size_t bin_size = (size_t)atoi(argv[4]);
    std::cout << "Bin size: " << bin_size << std::endl;
    
    uint64_t min_counter = (uint64_t)atoi(argv[5]);
    uint64_t max_counter = (uint64_t)atoi(argv[6]);
    std::cout << "Taking kmers of range " << min_counter << " to " << max_counter << std::endl;
    
    /* Construct histogram */
    // hist0 is total kmers; hist[1] is total distinct kmers (F0); from 2, hist[i] represents number of kmers with freq=i-1
    CKMCFileInfo kmc_info;
    kmc_database.Info(kmc_info);
    uint32_t k = kmc_info.kmer_length;
    
    std::cout << "KMC database read. k=" << kmc_info.kmer_length << std::endl;
    
    CKmerAPI * kmer_object = new CKmerAPI(kmc_info.kmer_length);
    
    gzFile fp = gzopen(argv[2], "r");
    kseq_t *seq;
    seq = kseq_init(fp);
    int l;
    
    std::ofstream output;
    output.open(argv[3]);
    
    
    while((l = kseq_read(seq)) >= 0) {
        std::cerr << "Processing: " << seq->name.s << std::endl;
        std::string current_sequence(seq->seq.s);
        
        
        for(size_t i = 0; i < seq->seq.l; i += bin_size) {
            uint64_t current_bin_total = 0;
            uint64_t count_below_min = 0;
            uint64_t count_above_max = 0;
            uint64_t count_n_kmer = 0;
            
            size_t max_idx = i + bin_size;
            
            int32_t last_n_pos = -1;
            
            std::vector<uint64_t> current_frequencies;
            
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
                
                uint64_t mode = current_mode;
                uint64_t median;
                if(current_frequencies.size() & 1) median = current_frequencies[current_frequencies.size() / 2];
                else if(current_frequencies.size() > 0) median = (current_frequencies[current_frequencies.size() / 2] + current_frequencies[current_frequencies.size() / 2 - 1]) / 2;
                
                output << seq->name.s << "\t" << i << "\t" << max_idx << "\t" << current_bin_total << "\t" << current_frequencies.size() << "\t" << double_t(current_bin_total)/double_t(current_frequencies.size()) << "\t" << mode << "\t" << median << "\t" << count_n_kmer << "\t" << count_below_min << "\t" << count_above_max << "\n";
            } else {
                uint64_t num_bases = max_idx - i;
                output << seq->name.s << "\t" << i << "\t" << max_idx << "\t" << 0 << "\t" << 0 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << count_n_kmer << "\t" << count_below_min << "\t" << count_above_max << "\n";
            }
        }
    }
    
    output.close();
    
    return 0;
}
