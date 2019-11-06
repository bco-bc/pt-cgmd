#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/file.o \
	${OBJECTDIR}/src/param.o \
	${OBJECTDIR}/src/poisson-process.o \
	${OBJECTDIR}/src/telegraph-process.o \
	${OBJECTDIR}/src/util.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f5 \
	${TESTDIR}/TestFiles/f1 \
	${TESTDIR}/TestFiles/f6 \
	${TESTDIR}/TestFiles/f9 \
	${TESTDIR}/TestFiles/f8 \
	${TESTDIR}/TestFiles/f4 \
	${TESTDIR}/TestFiles/f7 \
	${TESTDIR}/TestFiles/f3 \
	${TESTDIR}/TestFiles/f2

# Test Object Files
TESTOBJECTFILES= \
	${TESTDIR}/tests/box-cube-test.o \
	${TESTDIR}/tests/cvector-test.o \
	${TESTDIR}/tests/matrix-map-test.o \
	${TESTDIR}/tests/mu-units-test.o \
	${TESTDIR}/tests/nint-test.o \
	${TESTDIR}/tests/seed-value-test.o \
	${TESTDIR}/tests/telegraph-process-test.o \
	${TESTDIR}/tests/value-cvector-test.o \
	${TESTDIR}/tests/value_type_test.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcpputil.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcpputil.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcpputil.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/src/file.o: src/file.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/file.o src/file.cpp

${OBJECTDIR}/src/param.o: src/param.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/param.o src/param.cpp

${OBJECTDIR}/src/poisson-process.o: src/poisson-process.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/poisson-process.o src/poisson-process.cpp

${OBJECTDIR}/src/telegraph-process.o: src/telegraph-process.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/telegraph-process.o src/telegraph-process.cpp

${OBJECTDIR}/src/util.o: src/util.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/util.o src/util.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-tests-subprojects .build-conf ${TESTFILES}
.build-tests-subprojects:

${TESTDIR}/TestFiles/f5: ${TESTDIR}/tests/box-cube-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f5 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/cvector-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f6: ${TESTDIR}/tests/matrix-map-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f6 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f9: ${TESTDIR}/tests/mu-units-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f9 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f8: ${TESTDIR}/tests/nint-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f8 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f4: ${TESTDIR}/tests/seed-value-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f4 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f7: ${TESTDIR}/tests/telegraph-process-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f7 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f3: ${TESTDIR}/tests/value-cvector-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f3 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f2: ${TESTDIR}/tests/value_type_test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f2 $^ ${LDLIBSOPTIONS}   


${TESTDIR}/tests/box-cube-test.o: tests/box-cube-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/box-cube-test.o tests/box-cube-test.cpp


${TESTDIR}/tests/cvector-test.o: tests/cvector-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/cvector-test.o tests/cvector-test.cpp


${TESTDIR}/tests/matrix-map-test.o: tests/matrix-map-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/matrix-map-test.o tests/matrix-map-test.cpp


${TESTDIR}/tests/mu-units-test.o: tests/mu-units-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/mu-units-test.o tests/mu-units-test.cpp


${TESTDIR}/tests/nint-test.o: tests/nint-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/nint-test.o tests/nint-test.cpp


${TESTDIR}/tests/seed-value-test.o: tests/seed-value-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/seed-value-test.o tests/seed-value-test.cpp


${TESTDIR}/tests/telegraph-process-test.o: tests/telegraph-process-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/telegraph-process-test.o tests/telegraph-process-test.cpp


${TESTDIR}/tests/value-cvector-test.o: tests/value-cvector-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/value-cvector-test.o tests/value-cvector-test.cpp


${TESTDIR}/tests/value_type_test.o: tests/value_type_test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iinclude -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/value_type_test.o tests/value_type_test.cpp


${OBJECTDIR}/src/file_nomain.o: ${OBJECTDIR}/src/file.o src/file.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/file.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iinclude -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/file_nomain.o src/file.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/file.o ${OBJECTDIR}/src/file_nomain.o;\
	fi

${OBJECTDIR}/src/param_nomain.o: ${OBJECTDIR}/src/param.o src/param.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/param.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iinclude -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/param_nomain.o src/param.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/param.o ${OBJECTDIR}/src/param_nomain.o;\
	fi

${OBJECTDIR}/src/poisson-process_nomain.o: ${OBJECTDIR}/src/poisson-process.o src/poisson-process.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/poisson-process.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iinclude -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/poisson-process_nomain.o src/poisson-process.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/poisson-process.o ${OBJECTDIR}/src/poisson-process_nomain.o;\
	fi

${OBJECTDIR}/src/telegraph-process_nomain.o: ${OBJECTDIR}/src/telegraph-process.o src/telegraph-process.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/telegraph-process.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iinclude -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/telegraph-process_nomain.o src/telegraph-process.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/telegraph-process.o ${OBJECTDIR}/src/telegraph-process_nomain.o;\
	fi

${OBJECTDIR}/src/util_nomain.o: ${OBJECTDIR}/src/util.o src/util.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/util.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iinclude -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/util_nomain.o src/util.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/util.o ${OBJECTDIR}/src/util_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f5 || true; \
	    ${TESTDIR}/TestFiles/f1 || true; \
	    ${TESTDIR}/TestFiles/f6 || true; \
	    ${TESTDIR}/TestFiles/f9 || true; \
	    ${TESTDIR}/TestFiles/f8 || true; \
	    ${TESTDIR}/TestFiles/f4 || true; \
	    ${TESTDIR}/TestFiles/f7 || true; \
	    ${TESTDIR}/TestFiles/f3 || true; \
	    ${TESTDIR}/TestFiles/f2 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
