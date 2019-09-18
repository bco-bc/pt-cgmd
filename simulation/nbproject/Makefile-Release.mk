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
	${OBJECTDIR}/src/cg-pol-water.o \
	${OBJECTDIR}/src/interactor.o \
	${OBJECTDIR}/src/langevin-velocity-verlet.o \
	${OBJECTDIR}/src/leap-frog.o \
	${OBJECTDIR}/src/lj-coulomb-forces.o \
	${OBJECTDIR}/src/model-factory.o \
	${OBJECTDIR}/src/no-bc.o \
	${OBJECTDIR}/src/pbc.o \
	${OBJECTDIR}/src/sfactory.o \
	${OBJECTDIR}/src/sim-data.o \
	${OBJECTDIR}/src/sim-model.o \
	${OBJECTDIR}/src/simulation.o \
	${OBJECTDIR}/src/velocity-verlet.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1 \
	${TESTDIR}/TestFiles/f3 \
	${TESTDIR}/TestFiles/f4 \
	${TESTDIR}/TestFiles/f2

# Test Object Files
TESTOBJECTFILES= \
	${TESTDIR}/tests/displacer-test.o \
	${TESTDIR}/tests/model-factory-test.o \
	${TESTDIR}/tests/pair-list-test.o \
	${TESTDIR}/tests/simulation-test.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-pthread
CXXFLAGS=-pthread

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-Wl,-rpath,'../cpputil/dist/Debug/GNU-Linux' -L../cpputil/dist/Debug/GNU-Linux -lcpputil -Wl,-rpath,'../particles/dist/Debug/GNU-Linux' -L../particles/dist/Debug/GNU-Linux -lparticles

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libsimulation.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libsimulation.${CND_DLIB_EXT}: ../cpputil/dist/Debug/GNU-Linux/libcpputil.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libsimulation.${CND_DLIB_EXT}: ../particles/dist/Debug/GNU-Linux/libparticles.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libsimulation.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libsimulation.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -lpthread -pthread -shared -fPIC

${OBJECTDIR}/src/cg-pol-water.o: src/cg-pol-water.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/cg-pol-water.o src/cg-pol-water.cpp

${OBJECTDIR}/src/interactor.o: src/interactor.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/interactor.o src/interactor.cpp

${OBJECTDIR}/src/langevin-velocity-verlet.o: src/langevin-velocity-verlet.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/langevin-velocity-verlet.o src/langevin-velocity-verlet.cpp

${OBJECTDIR}/src/leap-frog.o: src/leap-frog.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/leap-frog.o src/leap-frog.cpp

${OBJECTDIR}/src/lj-coulomb-forces.o: src/lj-coulomb-forces.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lj-coulomb-forces.o src/lj-coulomb-forces.cpp

${OBJECTDIR}/src/model-factory.o: src/model-factory.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/model-factory.o src/model-factory.cpp

${OBJECTDIR}/src/no-bc.o: src/no-bc.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/no-bc.o src/no-bc.cpp

${OBJECTDIR}/src/pbc.o: src/pbc.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pbc.o src/pbc.cpp

${OBJECTDIR}/src/sfactory.o: src/sfactory.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sfactory.o src/sfactory.cpp

${OBJECTDIR}/src/sim-data.o: src/sim-data.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sim-data.o src/sim-data.cpp

${OBJECTDIR}/src/sim-model.o: src/sim-model.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sim-model.o src/sim-model.cpp

${OBJECTDIR}/src/simulation.o: src/simulation.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/simulation.o src/simulation.cpp

${OBJECTDIR}/src/velocity-verlet.o: src/velocity-verlet.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/velocity-verlet.o src/velocity-verlet.cpp

# Subprojects
.build-subprojects:
	cd ../cpputil && ${MAKE}  -f Makefile CONF=Debug
	cd ../particles && ${MAKE}  -f Makefile CONF=Debug

# Build Test Targets
.build-tests-conf: .build-tests-subprojects .build-conf ${TESTFILES}
.build-tests-subprojects:

${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/displacer-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f3: ${TESTDIR}/tests/model-factory-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f3 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f4: ${TESTDIR}/tests/pair-list-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f4 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f2: ${TESTDIR}/tests/simulation-test.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f2 $^ ${LDLIBSOPTIONS}   


${TESTDIR}/tests/displacer-test.o: tests/displacer-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/displacer-test.o tests/displacer-test.cpp


${TESTDIR}/tests/model-factory-test.o: tests/model-factory-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/model-factory-test.o tests/model-factory-test.cpp


${TESTDIR}/tests/pair-list-test.o: tests/pair-list-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/pair-list-test.o tests/pair-list-test.cpp


${TESTDIR}/tests/simulation-test.o: tests/simulation-test.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/simulation-test.o tests/simulation-test.cpp


${OBJECTDIR}/src/cg-pol-water_nomain.o: ${OBJECTDIR}/src/cg-pol-water.o src/cg-pol-water.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/cg-pol-water.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/cg-pol-water_nomain.o src/cg-pol-water.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/cg-pol-water.o ${OBJECTDIR}/src/cg-pol-water_nomain.o;\
	fi

${OBJECTDIR}/src/interactor_nomain.o: ${OBJECTDIR}/src/interactor.o src/interactor.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/interactor.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/interactor_nomain.o src/interactor.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/interactor.o ${OBJECTDIR}/src/interactor_nomain.o;\
	fi

${OBJECTDIR}/src/langevin-velocity-verlet_nomain.o: ${OBJECTDIR}/src/langevin-velocity-verlet.o src/langevin-velocity-verlet.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/langevin-velocity-verlet.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/langevin-velocity-verlet_nomain.o src/langevin-velocity-verlet.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/langevin-velocity-verlet.o ${OBJECTDIR}/src/langevin-velocity-verlet_nomain.o;\
	fi

${OBJECTDIR}/src/leap-frog_nomain.o: ${OBJECTDIR}/src/leap-frog.o src/leap-frog.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/leap-frog.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/leap-frog_nomain.o src/leap-frog.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/leap-frog.o ${OBJECTDIR}/src/leap-frog_nomain.o;\
	fi

${OBJECTDIR}/src/lj-coulomb-forces_nomain.o: ${OBJECTDIR}/src/lj-coulomb-forces.o src/lj-coulomb-forces.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/lj-coulomb-forces.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lj-coulomb-forces_nomain.o src/lj-coulomb-forces.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/lj-coulomb-forces.o ${OBJECTDIR}/src/lj-coulomb-forces_nomain.o;\
	fi

${OBJECTDIR}/src/model-factory_nomain.o: ${OBJECTDIR}/src/model-factory.o src/model-factory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/model-factory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/model-factory_nomain.o src/model-factory.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/model-factory.o ${OBJECTDIR}/src/model-factory_nomain.o;\
	fi

${OBJECTDIR}/src/no-bc_nomain.o: ${OBJECTDIR}/src/no-bc.o src/no-bc.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/no-bc.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/no-bc_nomain.o src/no-bc.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/no-bc.o ${OBJECTDIR}/src/no-bc_nomain.o;\
	fi

${OBJECTDIR}/src/pbc_nomain.o: ${OBJECTDIR}/src/pbc.o src/pbc.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pbc.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pbc_nomain.o src/pbc.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pbc.o ${OBJECTDIR}/src/pbc_nomain.o;\
	fi

${OBJECTDIR}/src/sfactory_nomain.o: ${OBJECTDIR}/src/sfactory.o src/sfactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/sfactory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sfactory_nomain.o src/sfactory.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/sfactory.o ${OBJECTDIR}/src/sfactory_nomain.o;\
	fi

${OBJECTDIR}/src/sim-data_nomain.o: ${OBJECTDIR}/src/sim-data.o src/sim-data.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/sim-data.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sim-data_nomain.o src/sim-data.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/sim-data.o ${OBJECTDIR}/src/sim-data_nomain.o;\
	fi

${OBJECTDIR}/src/sim-model_nomain.o: ${OBJECTDIR}/src/sim-model.o src/sim-model.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/sim-model.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sim-model_nomain.o src/sim-model.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/sim-model.o ${OBJECTDIR}/src/sim-model_nomain.o;\
	fi

${OBJECTDIR}/src/simulation_nomain.o: ${OBJECTDIR}/src/simulation.o src/simulation.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/simulation.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/simulation_nomain.o src/simulation.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/simulation.o ${OBJECTDIR}/src/simulation_nomain.o;\
	fi

${OBJECTDIR}/src/velocity-verlet_nomain.o: ${OBJECTDIR}/src/velocity-verlet.o src/velocity-verlet.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/velocity-verlet.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -I../particles/include -std=c++14 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/velocity-verlet_nomain.o src/velocity-verlet.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/velocity-verlet.o ${OBJECTDIR}/src/velocity-verlet_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f1 || true; \
	    ${TESTDIR}/TestFiles/f3 || true; \
	    ${TESTDIR}/TestFiles/f4 || true; \
	    ${TESTDIR}/TestFiles/f2 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} -r ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcpputil.so ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.so
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libsimulation.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:
	cd ../cpputil && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../particles && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
