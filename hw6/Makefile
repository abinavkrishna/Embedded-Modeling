# Makefile: Canny Edge Detector model in SystemC

# --- settings ---

VIDEO	= EngPlaza
FRAMES	= $(VIDEO)[0-9][0-9][0-9]_edges.pgm

SYSTEMC_HOME = /opt/pkg/systemc-2.3.1

CC	= g++
CCOPT	= -g -Wall

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
	rm -f *.o
	rm -f canny
	rm -f $(FRAMES)

# --- compile the example ---

canny: canny.cpp
	$(CC) $(CCOPT) $< \
		-I$(SYSTEMC_HOME)/include \
		-L$(SYSTEMC_HOME)/lib-linux64 \
		-Xlinker -R -Xlinker $(SYSTEMC_HOME)/lib-linux64 \
		-lsystemc -o $@

# EOF

