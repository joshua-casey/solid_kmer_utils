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
./reference_kmer_counter <kmc database path> <sequence file to scan> <output tsv file> <bin size>
```

For example:
```
./reference_kmer_counter kmc_result.res chm13.draft_v1.1.fasta test_out.tsv 10000
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
| 9 | count_0 | Number of kmers in bin with frequency 0 |
| 10 | count_non_0 | Number of kmers in bin with frequency > 0 |
| 11 | average_excl_0 | Average frequency of kmers, excluding 0-frequency (i.e. total_freq / count_non_0) (if empty bin, i.e. all Ns: -1) |
| 12 | mode_excl_0 | Mode of the kmer frequncies, excluding 0-frequency (if empty bin, i.e. all Ns: -1) |
| 13 | median_excl_0 | Median of the kmer frequencies, excluding 0-frequency (if empty bin, i.e. all Ns: -1) |
| 14 | count_N | Number of kmers with "N" |
