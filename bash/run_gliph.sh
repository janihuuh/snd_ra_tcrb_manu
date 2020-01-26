
## Run GLIPH (Glanville et al., Nature 2017) to cluster TCRbs to potentially epitope-specific clusters
## Here, we disrecard the motif search and focus only on TCRbs that differ by Hamming distance = 1

gliph=/Users/$me/Dropbox/aplastic_anemia_tcr/applications/gliph/gliph/bin/gliph-group-discovery.pl

## All TCRbs
$gliph --kmer_mindepth=1000000000 --simdepth=1 --gccutoff=1 --discontinuous=0 --tcr $files/results/gliph/input_files/dr_gliph_tot.txt
$gliph --kmer_mindepth=1000000000 --simdepth=1 --gccutoff=1 --discontinuous=0 --tcr $files/results/gliph/input_files/sn_gliph_tot.txt
$gliph --kmer_mindepth=1000000000 --simdepth=1 --gccutoff=1 --discontinuous=0 --tcr $files/results/gliph/input_files/sp_gliph_tot.txt
$gliph --kmer_mindepth=1000000000 --simdepth=1 --gccutoff=1 --discontinuous=0 --tcr $files/results/gliph/input_files/hc_gliph_tot.txt

## All non-naive TCRbs: exclude singletons
$gliph --kmer_mindepth=1000000000 --simdepth=1 --gccutoff=1 --discontinuous=0 --tcr $files/results/gliph/input_files/dr_gliph_min1.txt
$gliph --kmer_mindepth=1000000000 --simdepth=1 --gccutoff=1 --discontinuous=0 --tcr $files/results/gliph/input_files/sn_gliph_min1.txt
$gliph --kmer_mindepth=1000000000 --simdepth=1 --gccutoff=1 --discontinuous=0 --tcr $files/results/gliph/input_files/sp_gliph_min1.txt
$gliph --kmer_mindepth=1000000000 --simdepth=1 --gccutoff=1 --discontinuous=0 --tcr $files/results/gliph/input_files/hc_gliph_min1.txt
$gliph --kmer_mindepth=1000000000 --simdepth=1 --gccutoff=1 --discontinuous=0 --tcr $files/results/gliph/input_files/tot_gliph_min1.txt
