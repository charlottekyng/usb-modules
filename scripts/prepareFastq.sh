#!/bin/bash
# rename fastq files and create samples.txt and possibly samples.split.txt
mkdir -p fastq
find rawdata -name '*.gz' | xargs ln -t fastq

for x in fastq/*.fastq.gz; do
    if [ `grep -o "_" <<< "$x" | wc -l` -gt 1 ]; then
	mv $x `echo $x | sed 's/-//' | sed 's/^\([^_]*\)_\([^_]*\)_L\([0-9]*\)_R\([12]\)_\([0-9]*\)/\1-\2_\3\5.\4/' | sed 's/_//g' | sed 's/-/_/g'`
    fi;
done


paste <(ls fastq/*.fastq.gz | sed 's:.*/::; s/[._].*//') <(ls fastq/*.fastq.gz | sed 's:.*/::; s/\..*//') \
| awk 'BEGIN { OFS = "\t" } $1 != $2 { print }' | uniq > samples.split.txt
ls fastq/*.fastq.gz | sed 's:.*/::; s/_.*//' | sort | uniq >>samples.txt
