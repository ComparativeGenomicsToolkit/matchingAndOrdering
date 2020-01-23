#Modify this variable to set the location of sonLib
sonLibRootPath?=${rootPath}/../sonLib
sonLibPath=${sonLibRootPath}/lib
#Use sonLib bin and lib dirs
binPath=${sonLibRootPath}/bin
libPath=${sonLibPath}
export PATH := ${binPath}:${PATH}

include  ${sonLibRootPath}/include.mk

basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a ${dblibs}
basicLibsDependencies = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a 
