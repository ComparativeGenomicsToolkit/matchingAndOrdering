rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c
testBin = tests/testBin

all : externalTools ${libPath}/matchingAndOrdering.a ${binPath}/matchingAndOrderingTests ${testBin}/referenceMedianProblemTest

${libPath}/matchingAndOrdering.a : ${libSources} ${libHeaders} ${basicLibsDependencies} externalTools
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc matchingAndOrdering.a *.o
	ranlib matchingAndOrdering.a 
	rm *.o
	mv matchingAndOrdering.a ${libPath}/
	cp ${libHeaders} ${libPath}/

${binPath}/matchingAndOrderingTests : ${libTests} ${libSources} ${libHeaders} ${libPath}/matchingAndOrdering.a ${basicLibsDependencies} externalTools
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/matchingAndOrderingTests ${libTests} ${libPath}/matchingAndOrdering.a ${basicLibs}

${testBin}/referenceMedianProblemTest : ${testBin}/referenceMedianProblemTest.c ${libSources} ${libHeaders} ${libPath}/matchingAndOrdering.a ${basicLibsDependencies} externalTools
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${testBin}/referenceMedianProblemTest ${testBin}/referenceMedianProblemTest.c ${libPath}/matchingAndOrdering.a ${basicLibs}

clean : 
	rm -f *.o
	rm -f ${libPath}/matchingAndOrdering.a ${binPath}/matchingAndOrderingTests ${testBin}/referenceMedianProblemTest
	
test : all
	python allTests.py
