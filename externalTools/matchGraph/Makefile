include ../../include.mk
progs = mwmatching.py matchGraph.py

all: ${progs:%=${binPath}/%}

${binPath}/%: %
	cp -f $< $@
	chmod 775 $@

clean:
	rm -f  ${progs:%=${binPath}/%}
