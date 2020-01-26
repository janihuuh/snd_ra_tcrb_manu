
## Init objects to vdjtools scripts

############################

me=$(whoami)
application_files=/Users/$me/Dropbox/aplastic_anemia_tcr/
files=/Users/$me/Dropbox/ra_sn_tcrb/

vdj=$application_files/applications/vdjtools-1.2.1/vdjtools-1.2.1.jar
vdjdb=$application_files/applications/vdjdb-1.1.5/vdjdb-1.1.5.jar
vdjdb_new=$application_files/data/selected_tcrb/databases/vdjdb_new

############################

cd $files/data/tcr_filtered/
full=$(ls -d "$PWD"/*);

cd $files/data/tcr_filtered_10k/
over_10k=$(ls -d "$PWD"/*);

cd $files/data/tcr_filtered_10k_downsampled/
dn_10k=$(ls -d "$PWD"/*);

cd $files/data/tcr_filtered_10k_downsampled/
dr=$(ls -d "$PWD"/DR*);
sn=$(ls -d "$PWD"/SN*);
sp=$(ls -d "$PWD"/SP*);
hc=$(ls -d "$PWD"/HC*);

cd $files
clear

##############################



## Calc diversity stats for each sample with downsampled to 10k reads to normalize
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 10000 $full $files/results/diversity/unsampled_total
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 10000 $over_10k $files/results/diversity/unsampled_over10k

## Downsample individual samples to to 10k
java -Xmx4G -jar $vdj DownSample --size 10000 $over_10k $files/data/tcr_filtered_10k_downsampled/

## Find public clones, ie clones that are found at least in 2 people
java -Xmx4G -jar $vdj JoinSamples --intersect-type aa $dn_10k $files/results/public/downsampled_10k_

java -Xmx4G -jar $vdj JoinSamples --intersect-type aa $dr $hc $files/results/public/dr_hc
java -Xmx4G -jar $vdj JoinSamples --intersect-type aa $sn $hc $files/results/public/sn_hc
java -Xmx4G -jar $vdj JoinSamples --intersect-type aa $sp $hc $files/results/public/sp_hc

java -Xmx4G -jar $vdj JoinSamples --intersect-type aa $dr $sn $files/results/public/dr_sn
java -Xmx4G -jar $vdj JoinSamples --intersect-type aa $dr $sp $files/results/public/dr_sp
java -Xmx4G -jar $vdj JoinSamples --intersect-type aa $sp $sn $files/results/public/sp_sn

## Run vdjdb to find exact matches to previously known epitope-specific TCRs
java -Xmx4G -jar $vdjdb -R TRB -S human --vdjdb-conf 0 $files/results/tcrgp/input_files/dr_tcrgp_tot.txt $files/results/vdjdb/dr_
java -Xmx4G -jar $vdjdb -R TRB -S human --vdjdb-conf 0 $files/results/tcrgp/input_files/sn_tcrgp_tot.txt $files/results/vdjdb/sn_
java -Xmx4G -jar $vdjdb -R TRB -S human --vdjdb-conf 0 $files/results/tcrgp/input_files/sp_tcrgp_tot.txt $files/results/vdjdb/sp_
java -Xmx4G -jar $vdjdb -R TRB -S human --vdjdb-conf 0 $files/results/tcrgp/input_files/hc_tcrgp_tot.txt $files/results/vdjdb/hc_
