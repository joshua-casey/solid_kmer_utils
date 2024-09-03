#include "kmc_file.h"
#include "kmer_api.h"
#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
KSEQ_INIT(gzFile, gzread)

int main(int argc, char **argv) {   
    if(argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <KMC database path> <Reference path> <Output file> <Bin size>" << std::endl;
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
            uint64_t count_0 = 0;
            
            size_t max_idx = i + bin_size;
            if(max_idx >= seq->seq.l) max_idx = seq->seq.l;
            for(size_t j = i; j < max_idx; j++) {
                kmer_object->from_string_impl(current_sequence.begin() + j, k);
                uint64_t counter;
                if(kmc_database.CheckKmer(*kmer_object, counter)) {
                    current_bin_total += counter;
                } else {
                    count_0++;
                }
            }
            
            uint64_t num_bases = max_idx - i;
            output << seq->name.s << "\t" << i << "\t" << max_idx << "\t" << current_bin_total << "\t" << num_bases << "\t" << double_t(current_bin_total)/double_t(num_bases) << "\t" << count_0 << "\n";
        }
    }
    
    output.close();
    
    return 0;
}
