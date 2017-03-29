
##Set environment variables
#export TMPDIR = /scicore/scratch
#export PATH = $(HOME)/usr_nobackup:$(PATH)
#export R_LIBS := $(HOME)/R/x86_64-pc-linux-gnu-library/3.3

#LD_LIBRARY_PATH := $(HOME)/share/usr/mysql/lib:$(LD_LIBRARY_PATH)
#export LD_LIBRARY_PATH

#PERL5LIB := $(HOME)/share/usr/perl5/lib/perl5:$(HOME)/share/usr/ActivePerl-5.18/lib
#PERL5LIB := $(HOME)/share/usr/src/bioperl-live:$(PERL5LIB)
#PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl/modules
#PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-compara/modules
#PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-variation/modules
#PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-funcgen/modules
#PERL5LIB := $(HOME)/perl5/lib
#export PERL5LIB

#ifneq ("$(wildcard config.inc)", "")
#include config.inc
#endif
#ifneq ("$(wildcard project_config.inc)", "")
#include project_config.inc
#endif

# default reference:
# possible values: mm10, hg19, b37
#REF ?= b37
# possible values: ILLUMINA, IONTORRENT
#SEQ_PLATFORM ?= ILLUMINA
# possible values: NONE, BAITS, PCR
#CAPTURE_METHOD ?= NONE
# possible values: NONE, EXOME, IMPACTv5, CCP, CHPv2, IMPACT410, IMPACT341
#PANEL ?= NONE

SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt # generate using scripts/prepareFastq.sh [remove underscores from sample names]

#include usb-modules/genome_inc/$(REF).inc
#include usb-modules/config.inc

ifndef MAKEFILE_INC

PROJECT_NAME = $(shell pwd | sed 's:.*/projects/::; s:.*/data/::; s:.*kinglab/::; s:/:_:g')

ifneq ($(wildcard $(SAMPLE_FILE)),)
  SAMPLES ?= $(shell sed '/^\#/d' $(SAMPLE_FILE))
endif

get_tumors = $(wordlist 1,$(shell expr $(words $1) - 1),$1)
get_normal = $(lastword $1)
get_space_sets = $(shell sed '/^\#/d; s/\s\+/ /g; s/\s\+$$//;' $(SAMPLE_SET_FILE) | sed -n '$(1)p')
get_underscore_sets = $(shell sed '/^\#/d; s/\s\+/_/g; s/\s\+$$//;' $(SAMPLE_SET_FILE) | sed -n '$(1)p')

ifneq ($(wildcard $(SAMPLE_SET_FILE)),)
  NUM_SETS := $(shell sed '/^\#/d' $(SAMPLE_SET_FILE) | wc -l | cut -f 1 -d' ')
  SETS_SEQ := $(shell seq 1 $(NUM_SETS))
  $(foreach i,$(SETS_SEQ), \
	  $(eval set.$i := $(call get_space_sets,$i)))
  # define set_lookup.$(sample)
  $(foreach i,$(SETS_SEQ), \
	  $(foreach sample,$(set.$i), \
		  $(eval set_lookup.$(sample) := $i)))
  # define SAMPLE_SETS to contain sample_sample_normal
  $(foreach i,$(SETS_SEQ), \
	  $(eval SAMPLE_SETS += $(call get_underscore_sets,$i)))
  # tumor.SET => tumors
  $(foreach i,$(SETS_SEQ), \
	$(eval tumor.$(call get_underscore_sets,$i) := $(call get_tumors,$(set.$i))))
  # normal.SET => normal
  $(foreach i,$(SETS_SEQ), \
	$(eval normal.$(call get_underscore_sets,$i) := $(call get_normal,$(set.$i))))
  NORMAL_SAMPLES = $(filter-out poolednorm,$(foreach i,$(SETS_SEQ),$(call get_normal,$(set.$i))))
  TUMOR_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_tumors,$(set.$i)))
  SAMPLE_PAIRS = $(foreach set,$(SAMPLE_SETS),$(foreach tumor,$(tumor.$(set)),$(tumor)_$(normal.$(set))))
  # define tumor.pair and normal.pair to retrieve tumor/normal from pairs
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor.$(tumor)_$(normal.$(set)) := $(tumor))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval normal.$(tumor)_$(normal.$(set)) := $(normal.$(set)))))
  # tumor_normal.$(tumor) to retrieve tumor_normal
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor_normal.$(tumor) := $(tumor)_$(normal.$(set)))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor_normal.$(normal.$(set)) := $(tumor)_$(normal.$(set)))))
  # set.$(sample) to retrieve sample sets (underscores)
  $(foreach i,$(SETS_SEQ), \
	  $(foreach sample,$(set.$i), \
		  $(eval set.$(sample) := $(word $i,$(SAMPLE_SETS)))))
  UNMATCHED_SAMPLES = $(shell sed '/^\#/d' $(SAMPLE_FILE) $(SAMPLE_SET_FILE) | tr ' ' '\n' | sort | uniq -u)

  SAMPLE_SET_PAIRS = $(shell echo "$(SAMPLE_PAIRS) $(SAMPLE_SETS)" | tr ' ' '\n' | sort | uniq)
  $(foreach set,$(SAMPLE_SET_PAIRS), \
	$(eval samples.$(set) := $(shell echo "$(set)" | sed 's/_/ /g')))
endif

ifneq ($(wildcard $(SAMPLE_SPLIT_FILE)),)
  A = $(shell cut -f1 $(SAMPLE_SPLIT_FILE))
  B = $(shell cut -f2 $(SAMPLE_SPLIT_FILE))
  $(foreach i,$(shell seq 1 $(words $(A))),$(eval split.$(word $i,$(A)) += $(word $i,$(B))))
  UNSPLIT_SAMPLES = $(B)
  SPLIT_SAMPLES = $(shell cut -f1 $(SAMPLE_SPLIT_FILE) | sort | uniq)
endif

ALL_SAMPLES = $(SAMPLE_PAIRS) $(SAMPLES)

COUNT_SAMPLES = $(shell expr `sed 's/\..*//; s:.*/::' <<< $1 | grep -o "_" | wc -l` + 1)

NO_RM ?= false
ifeq ($(NO_RM),true)
    RM := touch
    RMR = touch
else
    RM := rm -f
    RMR := rm -r
endif

# define $(,) and $( ) for escaping commas and spaces
, := ,
space :=
space +=
$(space) := 
$(space) +=

# $(call strip-suffix,filename)
strip-suffix = $(firstword $(subst ., ,$(1)))

LINK = ln -svf $(notdir $1) $(2) && touch $1

NOW := $(shell date +"%F")

MKDIR = mkdir -p -m 775
MKDIRS = $(MKDIR) $(LOGDIR)/$(@D) $(@D)
LOG = $(LOGDIR)/$(@).log

# SGE variables
HOSTNAME = "*"
SGE_RREQ = -N X$(@F) -cwd -l hostname=$(HOSTNAME)
PARALLEL_ENV = smp
MEM_BY_CORE = -l membycore=$1
RUN_TIME = -l runtime=$1
REQ = SGE_RREQ="$(SGE_RREQ)"
REQ_PARALLEL_MEM = SGE_RREQ="$(SGE_RREQ) $(call MEM_BY_CORE,$2) -pe $(PARALLEL_ENV) $(1)"
REQ_MEM = SGE_RREQ="$(SGE_RREQ) $(call MEM_BY_CORE,$1) $(call RUN_TIME,$1)"
UMASK = 002
#ACTIVATE_JRFLAB_MODULES_ENV = unset PYTHONPATH; source $(JRFLAB_MODULES_ENV)/bin/activate $(JRFLAB_MODULES_ENV)
INIT = module purge; $(MKDIRS); umask $(UMASK); set -o pipefail; ####echo $(ACTIVATE_JRFLAB_MODULES_ENV);
INIT_NOCLUSTER = echo "module purge; umask $(UMASK); set -o pipefail"; #####$(ACTIVATE_JRFLAB_MODULES_ENV);";
INIT_MEM = module purge; $(MKDIRS); umask $(UMASK); $(call REQ_MEM,$1,$2);
INIT_PARALLEL_MEM = module purge; $(MKDIRS); umask $(UMASK); $(call REQ_PARALLEL_MEM,$1,$2,$3);

# openMPI
#MPIRUN = /opt/bin/mpirun
MPIRUN = /scicore/soft/apps/OpenMPI/1.8.4-GCC-4.8.4/bin/mpirun
MPI_ENV = openmpi
REQ_MPI_MEM = SGE_RREQ="$(SGE_RREQ) $(call MEM_BY_CORE,$2) -pe $(MPI_ENV) $(1)"
INIT_MPI_MEM = $(MKDIRS); umask $(UMASK); $(call REQ_MPI_MEM,$1,$2,$3); $(MPIRUN) -np $1

#QRSH = /common/sge/bin/lx-amd64/qrsh
QRSH = /scicore/soft/UGE/current/bin/lx-amd64/qrsh
#LAUNCH = $(INIT) $(QRSH) -N X$(@F) -V -cwd -q $(QUEUE) -now n -notify
#LAUNCH_PARALLEL_MEM = $(INIT) $(LAUNCH) -pe $(PARALLEL_ENV) $1 $(call MEM_FREE,$2,$3)
#LAUNCH_MEM = $(INIT) $(LAUNCH) $(call MEM_FREE,$1,$2)
#LAUNCH = $(INIT) $(QRSH) -N X$(@F) -V -cwd -now n -notify
LAUNCH = $(INIT) $(QRSH) -N X$(@F) -cwd -now n -notify
LAUNCH_PARALLEL_MEM = $(INIT) $(LAUNCH) -pe $(PARALLEL_ENV) $1 $(call MEM_BY_CORE,$2) $(call RUN_TIME $3)
LAUNCH_MEM = $(INIT) $(LAUNCH) $(call MEM_BY_CORE,$1)

# launch read script from stdin
QSUB_PRIORITY ?= 0
#override QSUB_ARGS = -V -wd $(PWD) -now n -notify -o $(LOG) -j y -S $(SHELL) -p $(QSUB_PRIORITY) -b n
override QSUB_ARGS = -wd $(PWD) -now n -notify -o $(LOG) -j y -S $(SHELL) -p $(QSUB_PRIORITY) -b n
QMAKE_BIN = /scicore/soft/UGE/current/bin/lx-amd64/qmake
QMAKE_QUEUE = $(QMAKE_BIN) -q $1
QMAKE = $(call QSUB_QUEUE,$(QUEUE))

QSUB_BIN_OPTS = -o $@
ifeq ($(USE_QSUB_DAEMON),true)
QSUB_BIN = usb-modules/scripts/qsubClient.pl
else
QSUB_BIN = usb-modules/scripts/qsub.pl
endif
QSUB_QUEUE = $(QSUB_BIN) $(QSUB_BIN_OPTS) --# -q $1
QSUB = $(QSUB_BIN) $(QSUB_BIN_OPTS) --# -q $(QUEUE)
# check output file
QSUB_CHECK = $(QSUB_BIN) $(QSUB_BIN_OPTS) -c --# -q $(QUEUE)
QSUB_CHECK_QUEUE = $(QSUB_BIN) $(QSUB_BIN_OPTS) -c --# -q $1

USE_CLUSTER ?= true
ifeq ($(USE_CLUSTER),false)
LSCRIPT = $(MKDIRS); ($(INIT_NOCLUSTER) echo $1) | $(SHELL) >& $(LOG)
LSCRIPT_PARALLEL_MEM = $(MKDIRS); ($(INIT_NOCLUSTER) echo $4) | $(SHELL) >& $(LOG)
LSCRIPT_MEM = $(MKDIRS); $ ($(INIT_NOCLUSTER) echo $3) | $(SHELL) >& $(LOG)

LSCRIPT_NET = $(LSCRIPT)
LSCRIPT_PARALLEL_MEM_NET = $(LSCRIPT_PARALLEL_MEM)
LSCRIPT_MEM_NET = $(LSCRIPT_MEM)

LSCRIPT_NAMED = $(call LSCRIPT,$2)
LSCRIPT_NAMED_PARALLEL_MEM = $(call LSCRIPT_PARALLEL_MEM,$2,$3,$4,$5)
LSCRIPT_NAMED_MEM = $(call LSCRIPT_MEM,$2,$3,$4)

LSCRIPT_CHECK = $(LSCRIPT)
LSCRIPT_CHECK_PARALLEL_MEM = $(LSCRIPT_PARALLEL_MEM)
LSCRIPT_CHECK_MEM = $(LSCRIPT_MEM)

LSCRIPT_CHECK_NET = $(LSCRIPT)
LSCRIPT_CHECK_PARALLEL_MEM_NET = $(LSCRIPT_PARALLEL_MEM)
LSCRIPT_CHECK_MEM_NET = $(LSCRIPT_MEM)

LSCRIPT_CHECK_NAMED = $(LSCRIPT_NAMED)
LSCRIPT_CHECK_NAMED_PARALLEL_MEM = $(LSCRIPT_NAMED_PARALLEL_MEM)
LSCRIPT_CHECK_NAMED_MEM = $(LSCRIPT_MEM)

else ## use cluster
LSCRIPT_NAMED = ($(INIT) echo "umask $(UMASK);" $2) | $(QSUB) -N $(PROJECT_NAME)_$1 $(QSUB_ARGS) $(call MEM_BY_CORE,1G) $(call RUN_TIME,00:29:59)
LSCRIPT_NAMED_PARALLEL_MEM = ($(INIT) echo "umask $(UMASK);" $5) | $(QSUB) -N $(PROJECT_NAME)_$1 $(QSUB_ARGS) -pe $(PARALLEL_ENV) $2 $(call MEM_BY_CORE,$3) $(call RUN_TIME,$4)
LSCRIPT_NAMED_MEM = ($(INIT) echo "umask $(UMASK);" $4) | $(QSUB) -N $(PROJECT_NAME)_$1 $(QSUB_ARGS) $(call MEM_BY_CORE,$2) $(call RUN_TIME,$3)

LSCRIPT = $(call LSCRIPT_NAMED_MEM,$(@F),1G,2G,$1)
LSCRIPT_PARALLEL_MEM = $(call LSCRIPT_NAMED_PARALLEL_MEM,$(@F),$1,$2,$3,$4)
LSCRIPT_MEM = $(call LSCRIPT_NAMED_MEM,$(@F),$1,$2,$3)

LSCRIPT_CHECK_NAMED = $(INIT) echo "umask $(UMASK);" $2 | $(QSUB_CHECK) -N $(PROJECT_NAME)_$1 $(QSUB_ARGS)
LSCRIPT_CHECK_NAMED_PARALLEL_MEM = ($(INIT) echo "umask $(UMASK);" $5) | $(QSUB_CHECK) -N $(PROJECT_NAME)_$1 $(QSUB_ARGS) -pe $(PARALLEL_ENV) $2 $(call MEM_BY_CORE,$3) $(call RUN_TIME,$4)
LSCRIPT_CHECK_NAMED_MEM = ($(INIT) echo "umask $(UMASK);" $4) | $(QSUB_CHECK) -N $(PROJECT_NAME)_$1 $(QSUB_ARGS) $(call MEM_BY_CORE,$2) $(call RUN_TIME,$3)

LSCRIPT_CHECK = $(call LSCRIPT_CHECK_NAMED,$(@F),$1)
LSCRIPT_CHECK_PARALLEL_MEM = $(call LSCRIPT_CHECK_NAMED_PARALLEL_MEM,$(@F),$1,$2,$3,$4)
LSCRIPT_CHECK_MEM = $(call LSCRIPT_CHECK_NAMED_MEM,$(@F),$1,$2,$3)

#LSCRIPT_NET = ($(INIT) echo "umask $(UMASK);" $1) | $(call QSUB_QUEUE,out.q) -N $(PROJECT_NAME)_$(@F) $(QSUB_ARGS) -l internet=1
#LSCRIPT_BIG = ($(INIT) echo "umask $(UMASK);" $1) | $(call QSUB_QUEUE,big.q) -N $(PROJECT_NAME)_$(@F) $(QSUB_ARGS)
##LSCRIPT_PARALLEL_MEM_NET = ($(INIT) echo "umask $(UMASK);" $4) | $(call QSUB_QUEUE,out.q) -N $(PROJECT_NAME)_$(@F) $(QSUB_ARGS) -pe $(PARALLEL_ENV) $1 $(call MEM_BY_CORE,$2) $(call RUN_TIME,$3) -l internet=1
#LSCRIPT_PARALLEL_MEM_BIG = ($(INIT) echo "umask $(UMASK);" $4) | $(call QSUB_QUEUE,big.q) -N $(PROJECT_NAME)_$(@F) $(QSUB_ARGS) -pe $(PARALLEL_ENV) $1 $(call MEM_BY_CORE,$2) $(call RUN_TIME,$3)
#LSCRIPT_MEM_NET = ($(INIT) echo "umask $(UMASK);" $3) | $(call QSUB_QUEUE,out.q) -N $(PROJECT_NAME)_$(@F) $(QSUB_ARGS) $(call MEM_BY_CORE,$1) $(call RUN_TIME,$2) -l internet=1
#LSCRIPT_MEM_BIG = ($(INIT) echo "umask $(UMASK);" $3) | $(call QSUB_QUEUE,big.q) -N $(PROJECT_NAME)_$(@F) $(QSUB_ARGS) $(call MEM_BY_CORE,$1) $(call RUNTIME,$2)
#LSCRIPT_CHECK_NET = ($(INIT) echo "umask $(UMASK);" $1) | $(call QSUB_CHECK_QUEUE,out.q) -N $(PROJECT_NAME)_$(@F) $(QSUB_ARGS) -l internet=1
#LSCRIPT_CHECK_PARALLEL_MEM_NET = ($(INIT) echo "umask $(UMASK);" $4) | $(call QSUB_CHECK_QUEUE,out.q) -N $(PROJECT_NAME)_$(@F) $(QSUB_ARGS) -pe $(PARALLEL_ENV) $1 $(call MEM_FREE,$2,$3) -l internet=1
#LSCRIPT_CHECK_MEM_NET = $(INIT) echo "umask $(UMASK);" $3 | $(call QSUB_CHECK_QUEUE,out.q) -N $(PROJECT_NAME)_$(@F) $(QSUB_ARGS) $(call MEM_FREE,$1,$2) -l internet=1

endif

#MD5 = md5sum $(1:.md5=) > $1
MD5 = md5sum $(@:.md5=) > $@
#CHECK_MD5 = md5sum -c $(filter %.md5,$^); i=0; until md5sum --status -c $(filter %.md5,$^); do if (( $$i % 5  == 0 )); then date; md5sum -c $(filter %.md5,$^); fi; if (( $$i > 20 )); then echo giving up; exit 1; fi; sleep 5; let i=i+1 ; done;
CHECK_MD5 = for i in {0..20}; do if md5sum -c $(filter %.md5,$^); then break; fi; sleep 5; done;
<M = $(<:.md5=)
^M = $(^:.md5=)
@M = $(@:.md5=)
<<M = $(word 2,$(^M))
<<<M = $(word 3,$(^M))

# $(call CHECK_VCF,vcf.file,target,command)
# check for valid header, ie non-empty vcf file
CHECK_VCF = if [ `grep -v '^\#' $1 | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $1 | wc -l` -eq 1 ]; then cp $1 $2; else $3; fi
CHECK_MAF = if [ `grep -v '^\#' $1 | wc -l` -eq 1 ] && [ `grep '^Hugo_Symbol' $1 | wc -l` -eq 1 ]; then cp $1 $2; else $3; fi
CHECK_VCF_CMD = if [ `grep -v '^\#' $1 | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $1 | wc -l` -eq 1 ]; then $2; else $3; fi

<< = $(word 2,$^)
<<< = $(word 3,$^)
<<<< = $(word 4,$^)
4< = $(word 4,$^)
<<<<< = $(word 5,$^)
5< = $(word 5,$^)

KNIT = $(RSCRIPT) modules/scripts/knit.R


define R_INIT
dir.create('$(@D)', showWarnings = F, recursive = T)
qw <- function(x) unlist(strsplit(x, "[[:space:]]+"))
endef

endif
MAKEFILE_INC = true
