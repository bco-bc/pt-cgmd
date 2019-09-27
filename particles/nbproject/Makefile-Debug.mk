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
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/atom.o \
	${OBJECTDIR}/src/atomistic.o \
	${OBJECTDIR}/src/bead.o \
	${OBJECTDIR}/src/coarse-grained.o \
	${OBJECTDIR}/src/continuous-protonatable-bead.o \
	${OBJECTDIR}/src/discrete-protonatable-bead.o \
	${OBJECTDIR}/src/particle-spec-catalog.o \
	${OBJECTDIR}/src/particle-spec.o \
	${OBJECTDIR}/src/particle.o \
	${OBJECTDIR}/src/pfactory.o \
	${OBJECTDIR}/src/protonation-site-catalog.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1 \
	${TESTDIR}/TestFiles/f2 \
	${TESTDIR}/TestFiles/f3 \
	${TESTDIR}/TestFiles/f5 \
	${TESTDIR}/TestFiles/f4 \
	${TESTDIR}/TestFiles/f6

# Test Object Files
TESTOBJECTFILES= \
	${TESTDIR}/tests/coarse-grained-test.o \
	${TESTDIR}/tests/particle-spec-catalog-test.o \
	${TESTDIR}/tests/particle-test.o \
	${TESTDIR}/tests/prot-site-catalog-test.o \
	${TESTDIR}/tests/protonation-site-test.o \
	${TESTDIR}/tests/read-cg-test.o

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
LDLIBSOPTIONS=-Wl,-rpath,'../cpputil/dist/Debug/GNU-Linux' -L../cpputil/dist/Debug/GNU-Linux -lcpputil

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.${CND_DLIB_EXT}: ../cpputil/dist/Debug/GNU-Linux/libcpputil.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/src/atom.o: src/atom.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/atom.o src/atom.cpp

${OBJECTDIR}/src/atomistic.o: src/atomistic.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/atomistic.o src/atomistic.cpp

${OBJECTDIR}/src/bead.o: src/bead.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bead.o src/bead.cpp

${OBJECTDIR}/src/coarse-grained.o: src/coarse-grained.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/coarse-grained.o src/coarse-grained.cpp

${OBJECTDIR}/src/continuous-protonatable-bead.o: src/continuous-protonatable-bead.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/continuous-protonatable-bead.o src/continuous-protonatable-bead.cpp

${OBJECTDIR}/src/discrete-protonatable-bead.o: src/discrete-protonatable-bead.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/discrete-protonatable-bead.o src/discrete-protonatable-bead.cpp

${OBJECTDIR}/src/particle-spec-catalog.o: src/particle-spec-catalog.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/particle-spec-catalog.o src/particle-spec-catalog.cpp

${OBJECTDIR}/src/particle-spec.o: src/particle-spec.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/particle-spec.o src/particle-spec.cpp

${OBJECTDIR}/src/particle.o: src/particle.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/particle.o src/particle.cpp

${OBJECTDIR}/src/pfactory.o: src/pfactory.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pfactory.o src/pfactory.cpp

${OBJECTDIR}/src/protonation-site-catalog.o: src/protonation-site-catalog.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/protonation-site-catalog.o src/protonation-site-catalog.cpp

# Subprojects
.build-subprojects:
	cd ../cpputil && ${MAKE}  -f Makefile CONF=Debug

# Build Test Targets
.build-tests-conf: .build-tests-subprojects .build-conf ${TESTFILES}
.build-tests-subprojects:

${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/coarse-grained-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f2: ${TESTDIR}/tests/particle-spec-catalog-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f2 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f3: ${TESTDIR}/tests/particle-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f3 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f5: ${TESTDIR}/tests/prot-site-catalog-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f5 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f4: ${TESTDIR}/tests/protonation-site-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f4 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f6: ${TESTDIR}/tests/read-cg-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f6 $^ ${LDLIBSOPTIONS}   


${TESTDIR}/tests/coarse-grained-test.o: tests/coarse-grained-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/coarse-grained-test.o tests/coarse-grained-test.cpp


${TESTDIR}/tests/particle-spec-catalog-test.o: tests/particle-spec-catalog-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/particle-spec-catalog-test.o tests/particle-spec-catalog-test.cpp


${TESTDIR}/tests/particle-test.o: tests/particle-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/particle-test.o tests/particle-test.cpp


${TESTDIR}/tests/prot-site-catalog-test.o: tests/prot-site-catalog-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/prot-site-catalog-test.o tests/prot-site-catalog-test.cpp


${TESTDIR}/tests/protonation-site-test.o: tests/protonation-site-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/protonation-site-test.o tests/protonation-site-test.cpp


${TESTDIR}/tests/read-cg-test.o: tests/read-cg-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/read-cg-test.o tests/read-cg-test.cpp


${OBJECTDIR}/src/atom_nomain.o: ${OBJECTDIR}/src/atom.o src/atom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/atom.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/atom_nomain.o src/atom.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/atom.o ${OBJECTDIR}/src/atom_nomain.o;\
	fi

${OBJECTDIR}/src/atomistic_nomain.o: ${OBJECTDIR}/src/atomistic.o src/atomistic.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/atomistic.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/atomistic_nomain.o src/atomistic.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/atomistic.o ${OBJECTDIR}/src/atomistic_nomain.o;\
	fi

${OBJECTDIR}/src/bead_nomain.o: ${OBJECTDIR}/src/bead.o src/bead.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/bead.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bead_nomain.o src/bead.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/bead.o ${OBJECTDIR}/src/bead_nomain.o;\
	fi

${OBJECTDIR}/src/coarse-grained_nomain.o: ${OBJECTDIR}/src/coarse-grained.o src/coarse-grained.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/coarse-grained.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/coarse-grained_nomain.o src/coarse-grained.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/coarse-grained.o ${OBJECTDIR}/src/coarse-grained_nomain.o;\
	fi

${OBJECTDIR}/src/continuous-protonatable-bead_nomain.o: ${OBJECTDIR}/src/continuous-protonatable-bead.o src/continuous-protonatable-bead.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/continuous-protonatable-bead.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/continuous-protonatable-bead_nomain.o src/continuous-protonatable-bead.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/continuous-protonatable-bead.o ${OBJECTDIR}/src/continuous-protonatable-bead_nomain.o;\
	fi

${OBJECTDIR}/src/discrete-protonatable-bead_nomain.o: ${OBJECTDIR}/src/discrete-protonatable-bead.o src/discrete-protonatable-bead.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/discrete-protonatable-bead.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/discrete-protonatable-bead_nomain.o src/discrete-protonatable-bead.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/discrete-protonatable-bead.o ${OBJECTDIR}/src/discrete-protonatable-bead_nomain.o;\
	fi

${OBJECTDIR}/src/particle-spec-catalog_nomain.o: ${OBJECTDIR}/src/particle-spec-catalog.o src/particle-spec-catalog.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/particle-spec-catalog.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/particle-spec-catalog_nomain.o src/particle-spec-catalog.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/particle-spec-catalog.o ${OBJECTDIR}/src/particle-spec-catalog_nomain.o;\
	fi

${OBJECTDIR}/src/particle-spec_nomain.o: ${OBJECTDIR}/src/particle-spec.o src/particle-spec.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/particle-spec.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/particle-spec_nomain.o src/particle-spec.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/particle-spec.o ${OBJECTDIR}/src/particle-spec_nomain.o;\
	fi

${OBJECTDIR}/src/particle_nomain.o: ${OBJECTDIR}/src/particle.o src/particle.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/particle.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/particle_nomain.o src/particle.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/particle.o ${OBJECTDIR}/src/particle_nomain.o;\
	fi

${OBJECTDIR}/src/pfactory_nomain.o: ${OBJECTDIR}/src/pfactory.o src/pfactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pfactory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pfactory_nomain.o src/pfactory.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pfactory.o ${OBJECTDIR}/src/pfactory_nomain.o;\
	fi

${OBJECTDIR}/src/protonation-site-catalog_nomain.o: ${OBJECTDIR}/src/protonation-site-catalog.o src/protonation-site-catalog.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/protonation-site-catalog.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/protonation-site-catalog_nomain.o src/protonation-site-catalog.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/protonation-site-catalog.o ${OBJECTDIR}/src/protonation-site-catalog_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f1 || true; \
	    ${TESTDIR}/TestFiles/f2 || true; \
	    ${TESTDIR}/TestFiles/f3 || true; \
	    ${TESTDIR}/TestFiles/f5 || true; \
	    ${TESTDIR}/TestFiles/f4 || true; \
	    ${TESTDIR}/TestFiles/f6 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} -r ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcpputil.so
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:
	cd ../cpputil && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
