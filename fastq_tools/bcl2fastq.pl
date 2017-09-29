run=$1
module load bcl2fastq/2.17.1.14-goolf-1.7.20
bcl2fastq -R intensities/${run} -o intensities/${run}/Data/Fastq/ -d 4 -p 4 --barcode-mismatches 0
