# This module runs ContEst on snp vcf files from gatk
# Author: inodb

##### MAKE INCLUDES #####
include modules/Makefile.inc

LOGDIR = log/contest.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: contest

contest : contest/snp_vcf/all_contamination.txt $(foreach sample,$(SAMPLES),contest/snp_vcf/$(sample)_contamination.txt)

# ContEst on gatk snp_vcf folder
contest/snp_vcf/%_contamination.txt : bam/%.bam snp_vcf/%.snps.vcf
	$(call LSCRIPT_MEM,7G,8G,"$(MKDIR) contest/snp_vcf; $(call CONTEST_MEM,2G) -I $(<) -B:genotypes$(,)vcf $(<<) -o $(@)")

contest/snp_vcf/all_contamination.txt : $(foreach sample,$(SAMPLES),contest/snp_vcf/$(sample)_contamination.txt)
	( \
		head -1 $< | sed "s/^/sample\t/"; \
		for s in $(^); do \
			grep -P "META\t" $$s | sed "s/^/`basename $$s _contamination.txt`/"; \
		done | sort -rnk5,5; \
	) > $@
