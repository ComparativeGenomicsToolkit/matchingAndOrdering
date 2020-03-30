rootPath = ../..
include ${rootPath}/include.mk

DIRS := . MinCost GEOM 

SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

LIBS :=  -lrt 
INCLUDES := 
LIBDIR := 

# Mac OS X specific stuff
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
	CXXFLAGS += -DPM_TIMER_GETRUSAGE
	LIBS := 
endif

all: ${BINDIR}/blossom5 ${BINDIR}/blossom_randGraph.py

${BINDIR}/blossom_randGraph.py: blossom_randGraph.py
	cp blossom_randGraph.py ${BINDIR}
	chmod +x ${BINDIR}/blossom_randGraph.py

${BINDIR}/blossom5: ${OBJS}
	${CXX} ${CXXFLAGS} ${LIBDIR} -o $@ ${OBJS} ${LIBS}
	chmod +x ${BINDIR}/blossom5

.cpp.o:
	$(CXX) $(CXXFLAGS) ${INCLUDES} $< -c -o $@

clean:
	rm -f ${OBJS}
	rm -f ${BINDIR}/blossom5 ${BINDIR}/blossom_randGraph.py
