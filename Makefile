rootPath = .
include ${rootPath}/include.mk

export PYTHONPATH = ../sonLib/src:..

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c
testBin = tests/testBin

CPPFLAGS += -Iimpl -I../sonLib/externalTools/cutest/

all : externalToolsM ${LIBDIR}/matchingAndOrdering.a ${BINDIR}/matchingAndOrderingTests ${testBin}/referenceMedianProblemTest2

externalToolsM : 
	cd externalTools && ${MAKE} all

${LIBDIR}/matchingAndOrdering.a : ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -c ${libSources}
	${AR} rc matchingAndOrdering.a *.o
	${RANLIB} matchingAndOrdering.a 
	mv matchingAndOrdering.a ${LIBDIR}/
	cp ${libHeaders} ${LIBDIR}/

${BINDIR}/matchingAndOrderingTests : ${libTests} ${libSources} ${libHeaders} ${LIBDIR}/matchingAndOrdering.a ${LIBDEPENDS}
	${CC} ${CFLAGS} ${CPPFLAGS} ${LDFLAGS} -o ${BINDIR}/matchingAndOrderingTests ${libTests} ${LIBDIR}/matchingAndOrdering.a ${LDLIBS}

${testBin}/referenceMedianProblemTest2 : ${testBin}/referenceMedianProblemTest2.c ${libSources} ${libHeaders} ${LIBDIR}/matchingAndOrdering.a ${LIBDEPENDS} 
	${CC} ${CFLAGS} ${CPPFLAGS} ${LDFLAGS} -o ${testBin}/referenceMedianProblemTest2 ${testBin}/referenceMedianProblemTest2.c ${LIBDIR}/matchingAndOrdering.a ${LDLIBS}

clean : 
	cd externalTools && ${MAKE} clean
	rm -f *.o
	rm -f ${LIBDIR}/matchingAndOrdering.a ${BINDIR}/matchingAndOrderingTests ${testBin}/referenceMedianProblemTest

test : all
	${PYTHON} allTests.py
