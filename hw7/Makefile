# Makefile: Canny Edge Detector model for profiling
# EECS 222, W20, A7 (combined with A6)
#
# (used cannyA6_ref.sc and 'scc -psi -pui -sl' to generated canny.cc, canny.h;
#  then combined/edited both to canny.cpp for profiling with GPROF

# --- settings ---

VIDEO	= EngPlaza
FRAMES	= $(VIDEO)[0-9][0-9][0-9]_edges.pgm

SPECC	= /opt/sce
SCC	= scc
SCCOPT	= -g -ww -vvv -psi -pui -sl

# --- targets ---

all:	canny

test:	canny
	./canny
	set -e;	\
		for f in video/$(FRAMES); do \
		diff `basename $$f` $$f; \
		done

clean:
	rm -f *~ *.bak *.BAK
	rm -f *.o *.h *.cc *.si *.sir
	rm -f canny
	rm -f $(FRAMES)

# --- compile the example ---

canny: canny.cpp
	g++ -g -pg -Wall -U__GNUC__ -isystem $(SPECC)/inc canny.cpp -L$(SPECC)/lib -Xlinker -R -Xlinker $(SPECC)/lib -lbit4 -lbit -lsim -llonglong -lqthread -o canny

# EOF
