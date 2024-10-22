include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/mutsigcv.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : mutsigcv

mutsigcv : mutsigcv/mutsigcv.sig_genes.txt
mutsigcv/mutsigcv.sig_genes.txt : mutsigcv/mutsigcv_input.maf
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_SHORT),"$(MUTSIGCV) $(MCR) $^ \
	$(MUTSIGCV_COVERAGE_REF) $(MUTSIGCV_COV_REF) mutsigcv $(MUTSIGCV_DICT_REF) $(MUTSIGCV_SEQ_REF_DIR) && \
	mv mutsigcv.* mutsigcv")
