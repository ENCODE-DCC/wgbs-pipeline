#==================================================================================================
# PROJECT: gemBS
# FILE: Makefile
# DATE: 01/11/2015
# AUTHOR(S): Marcos Fernandez <mfernandez@cnag.crg.eu>
# DESCRIPTION: Builds the external tools for gemBS
#==================================================================================================

# Definitions

CPLUS=g++

FOLDER_BIN=../bin

TOOLS=readNameClean

TOOLS_SRC=$(addsuffix .cpp, $(TOOLS))
TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(TOOLS))

LIBS:=-lm

ifeq ($(HAVE_ZLIB),1)
LIBS:=$(LIBS) -lz
endif
ifeq ($(HAVE_BZLIB),1)
LIBS:=$(LIBS) -lbz2
endif


# Flags
ARCH_FLAGS = -D__LINUX__
HAVE_ZLIB = 1
HAVE_BZLIB = 1
HAVE_OPENMP = 1

GENERAL_FLAGS=-fPIC -Wall
ifeq ($(HAVE_ZLIB),1)
GENERAL_FLAGS:=$(GENERAL_FLAGS) -DHAVE_ZLIB
endif
ifeq ($(HAVE_BZLIB),1)
GENERAL_FLAGS:=$(GENERAL_FLAGS) -DHAVE_BZLIB
endif
ifeq ($(HAVE_OPENMP),1)
GENERAL_FLAGS:=$(GENERAL_FLAGS) -DHAVE_OPENMP
endif

SUPPRESS_CHECKS=-DNDEBUG -DGT_NO_CONSISTENCY_CHECKS
DEBUG_FLAGS=-g -ggdb3 -rdynamic -DGT_DEBUG
OPTIMIZTION_FLAGS=-O4

all: TOOLS_FLAGS=-O4 $(GENERAL_FLAGS) $(SUPPRESS_CHECKS)
all: $(TOOLS_BIN)

static: TOOLS_FLAGS=-O4 $(GENERAL_FLAGS) $(SUPPRESS_CHECKS) -static
static: $(TOOLS_BIN)

debug: TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(DEBUG_FLAGS)
debug: $(TOOLS_BIN)

$(TOOLS_BIN): $(TOOLS_SRC)
	$(CPLUS) $(TOOLS_FLAGS) -o $@ $(notdir $@).cpp $(LIBS) 

