#==================================================================================================
# PROJECT: gemBS
# FILE: Makefile
# DATE: 01/11/2015
# AUTHOR(S): Marcos Fernandez <mfernandez@cnag.crg.eu>
# DESCRIPTION: Builds the external tools for gemBS
#==================================================================================================


# Definitions

CC=gcc

ROOT_PATH=..

TOOLS=filter_vcf
FOLDER_BIN=../bin/
TOOLS_SRC=$(addsuffix .c, $(TOOLS))
TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(TOOLS))
LOKI_LIBS:=-I../loki/include -L../loki/libsrc
LIBS:=-lgen -lm


ifeq ($(HAVE_ZLIB),1)
LIBS:=$(LIBS) -lz
endif
ifeq ($(HAVE_BZLIB),1)
LIBS:=$(LIBS) -lbz2
endif

all: TOOLS_FLAGS=-O3 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: $(TOOLS_BIN)
	
static: TOOLS_FLAGS=-O3 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS) -static
static: $(TOOLS_BIN)

debug: TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(TOOLS_BIN)

$(TOOLS_BIN): $(TOOLS_SRC)
	$(CC) --std=c99  $(TOOLS_FLAGS) -o $@ $(notdir $@).c $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LOKI_LIBS) $(LIBS) $(EXTRA_LIBS)
