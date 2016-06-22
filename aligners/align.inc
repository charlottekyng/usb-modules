ifndef ALIGN_INC

BAM_SUFFIX = $(subst $( ),.,$(strip \
        $(if $(findstring false,$(BAM_NO_SORT)),sorted)\
        $(if $(findstring false,$(BAM_NO_FILTER)),filtered)\
        $(if $(findstring true,$(BAM_FIX_RG)),rg)\
        $(if $(findstring rmdup,$(BAM_DUP_TYPE)),rmdup)\
        $(if $(findstring markdup,$(BAM_DUP_TYPE)),markdup)\
        $(if $(findstring false,$(BAM_NO_REALN)),realn)\
        $(if $(findstring false,$(BAM_NO_RECAL)),recal)\
        bam))
$(info BAM_SUFFIX-ALIGN is $(BAM_SUFFIX))

endif
ALIGN_INC = true
