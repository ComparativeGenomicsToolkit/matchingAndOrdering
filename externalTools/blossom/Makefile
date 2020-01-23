rootPath = ../..
include ${rootPath}/include.mk

DIRS := . MinCost GEOM 

SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

CXXFLAGS := ${cppflags}
CXX := ${cpp}
LIBS :=  -lrt 
INCLUDES := 
LIBDIR := 

# Mac OS X specific stuff
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
	CXXFLAGS += -DPM_TIMER_GETRUSAGE
	LIBS := 
endif

all: ${binPath}/blossom5 ${binPath}/blossom_randGraph.py

${binPath}/blossom_randGraph.py: blossom_randGraph.py
	cp blossom_randGraph.py ${binPath}
	chmod +x ${binPath}/blossom_randGraph.py

${binPath}/blossom5: ${OBJS}
	${CXX} ${CXXFLAGS} ${LIBDIR} -o $@ ${OBJS} ${LIBS}
	chmod +x ${binPath}/blossom5

.cpp.o:
	$(CXX) $(CXXFLAGS) ${INCLUDES} $< -c -o $@

clean:
	rm -f ${OBJS}
	rm -f ${binPath}/blossom5 ${binPath}/blossom_randGraph.py
