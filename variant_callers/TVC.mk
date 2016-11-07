include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/variantCaller.inc

LOGDIR ?= log/tvc.$(NOW)

VPATH ?= bam
VARIANT_TYPES ?= 

tvc_vcfs : $(foreach t
