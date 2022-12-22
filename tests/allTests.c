/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"

extern CuSuite* referenceProblem2TestSuite(void);
extern CuSuite* cyclesConstrainedMatchingAlgorithmsTestSuite(void);
extern CuSuite* matchingAlgorithmsTestSuite(void);

int referenceRunAllTests(void) {
    CuString *output = CuStringNew();
    CuSuite* suite = CuSuiteNew();
    CuSuiteAddSuite(suite, referenceProblem2TestSuite());
    CuSuiteAddSuite(suite, cyclesConstrainedMatchingAlgorithmsTestSuite());
    CuSuiteAddSuite(suite, matchingAlgorithmsTestSuite());

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    int i = suite->failCount > 0;
    CuSuiteDelete(suite);
    return i;
}

int main(int argc, char *argv[]) {
    if(argc == 2) {
        st_setLogLevelFromString(argv[1]);
    }
    int i = referenceRunAllTests();

    return i;
}
