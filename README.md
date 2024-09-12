Various solid kmers utility

## Compilation

```
./install_htslib.sh
mkdir build
cd build
cmake ..
make
```

### Getting average kmers frequencies from a list of sequence

First, run KMC:
```
kmc -k<k> -m<mem> @reads.txt kmc_result.res <tmp>
```
This will output two files `kmc_result.res.kmc_pre` and `kmc_result.res.kmc_suf`

Then, run the compiled program
```
./reference_kmer_counter -i <kmc database> -r <reference file> [optional arguments]
```
The list of arguments are as follows:
```
Mandatory args:
  -i, --kmc_db <str> 	 	 KMC database path (without pre/suf).
  -r, --reference <str> 	 The fasta file of the reference to scan (accepts gzipped file).
Optional args:
  -o, --output_prefix <str>  Output prefix: output will be PREFIXkmer_count.tsv and PREFIXkmer_histo.tsv [DEAFULT: ]
  -b, --bin_size <int> 	 	 Bin size [DEAFULT: 10000]
  -s, --bin_shift <int> 	 Bin sliding window shift [DEAFULT: 10000]
  -m, --min_kmer <int> 	 	 Min kmer count [DEAFULT: 2]
  -M, --max_kmer <int> 	 	 Max kmer count [DEAFULT: 200]
  -h, --help 	 	 	     Prints the usage.
```
Note that minimum, maximum kmer count, and the statistics will obey the KMC counters, i.e. set -ci, -cs, and -cx accordingly.

For example:
```
./reference_kmer_counter -i kmc_result.res -r chm13.draft_v1.1.fasta
```

The output `test_out.tsv` will be of the format:
| Column no  | Name | Description  |
|---|---|---|
| 1 | seq_name | Sequence name   |
| 2 | bin_start | Bin start (0 based, inclusive) |
| 3 | bin_end | Bin end (0 based, exclusive)  |
| 4 | total_freq | Total kmer frequency in bin |
| 5 | total_kmer | Total kmers in bin |
| 6 | average | Average frequency of kmers (i.e. total_freq / total_kmer). (if empty bin, i.e. all Ns: -1) |
| 7 | mode | Mode of the kmer frequncies (if empty bin, i.e. all Ns: -1) |
| 8 | median | Median of the kmer frequencies (if empty bin, i.e. all Ns: -1) |
| 9 | count_N | Number of kmers with "N" |
| 10 | below_min | Number of kmers with frequencies below the specified minimum frequency |
| 11 | above_max | Number of kmers with frequencies above the specified maximum frequency |

The output `hist_out.tsv` will be of the format:
| Column no  | Name | Description  |
|---|---|---|
| 1 | seq_name | Sequence name   |
| 2 | bin_start | Bin start (0 based, inclusive) |
| 3 | bin_end | Bin end (0 based, exclusive)  |
| 4+ | histogram value | Number of kmers with specific frequency, column no 4 will represent kmers with min_kmer_count, column no 5 kmers with min_kmer_count+1, etc. |
