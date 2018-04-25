SHELL = /bin/sh

LEX = flex -I
YACC = bison -y
YFLAGS = -d
CC = gcc
MY_CFLAGS = -g -O2 -W -Wall -pedantic -Wno-long-long -ffloat-store
LDFLAGS = 
LIBS = -lm 
CP = /bin/cp
MV = /bin/mv
CHMOD = /bin/chmod
prefix = /Volumes/Home/Users/heath/code/loki
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin
INSTALL_PERM = 755
EXTRA_CFLAGS = #-DDEBUG -DTRACE_PEEL
#DMALLOC_INC = -I/opt/local/include
#DMALLOC_LIB = -L/opt/local/lib]
#DMALLOC_FLAGS = -DUSE_DMALLOC $(DMALLOC_INC)

INCLUDES = -I../include -I../lokisrc/
CFLAGS = $(MY_CFLAGS) $(EXTRA_CFLAGS) $(INCLUDES) $(DMALLOC_FLAGS)
ALL_LIBS = $(DMALLOC_LIB) -L../lokisrc -L../libsrc -lloki -lgen $(LIBS)

FENRIS_SRC = alloc_fenris.c fenris_setup.c count_bits.c estimate_freq.c \
peel_freq.c fenris_pen.c pen_model.c peel_init.c fenris_utils.c \
fenris_simple_peel.c fenris_complex_peel.c calc_pen.c fseg_pen.c fenris.c

FENRIS_OBJ = ${SCAN_SRC:.c=.o} ${FENRIS_SRC:.c=.o}

PROGS = fenris

all: $(PROGS)

fenris: $(FENRIS_OBJ) ../libsrc/libgen.a ../include/config.h
	$(CC) $(CFLAGS) -o $@ $(FENRIS_OBJ) $(LDFLAGS) $(ALL_LIBS)
	
../bin:
	mkdir ../bin

install: $(PROGS) ../bin
	for file in $(PROGS); do \
	$(CP) $$file $(bindir)/; \
	$(CHMOD) $(INSTALL_PERM) $(bindir)/$$file ; \
	done

clean: 
	rm -f *~ *.o *.a *.bak a.out core seedfile depend

distclean: clean
	rm -f $(PROGS) Makefile

depend: $(FENRIS_SRC)
	/usr/X11R6/bin/makedepend $(INCLUDES) $(DMALLOC_FLAGS) $(FENRIS_SRC)
	touch depend
	
# DO NOT DELETE THIS LINE -- make depend depends on it.
