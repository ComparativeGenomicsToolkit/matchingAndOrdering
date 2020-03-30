rootPath = ../..
include ${rootPath}/include.mk
progs = mwmatching.py matchGraph.py

all: ${progs:%=${BINDIR}/%}

${BINDIR}/%: %
	cp -f $< $@
	chmod 775 $@

clean:
	rm -f  ${progs:%=${BINDIR}/%}
