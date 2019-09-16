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
	${OBJECTDIR}/src/atomistic.o \
	${OBJECTDIR}/src/bead.o \
	${OBJECTDIR}/src/coarse-grained.o \
	${OBJECTDIR}/src/particle-spec-catalog.o \
	${OBJECTDIR}/src/particle-spec.o \
	${OBJECTDIR}/src/particle.o \
	${OBJECTDIR}/src/pfactory.o \
	${OBJECTDIR}/src/protonatable-bead.o \
	${OBJECTDIR}/src/protonation-site-catalog.o


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
LDLIBSOPTIONS=-Wl,-rpath,'../cpputil/dist/Release/GNU-Linux' -L../cpputil/dist/Release/GNU-Linux -lcpputil

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.${CND_DLIB_EXT}: ../cpputil/dist/Release/GNU-Linux/libcpputil.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/src/atomistic.o: src/atomistic.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/atomistic.o src/atomistic.cpp

${OBJECTDIR}/src/bead.o: src/bead.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bead.o src/bead.cpp

${OBJECTDIR}/src/coarse-grained.o: src/coarse-grained.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/coarse-grained.o src/coarse-grained.cpp

${OBJECTDIR}/src/particle-spec-catalog.o: src/particle-spec-catalog.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/particle-spec-catalog.o src/particle-spec-catalog.cpp

${OBJECTDIR}/src/particle-spec.o: src/particle-spec.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/particle-spec.o src/particle-spec.cpp

${OBJECTDIR}/src/particle.o: src/particle.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/particle.o src/particle.cpp

${OBJECTDIR}/src/pfactory.o: src/pfactory.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pfactory.o src/pfactory.cpp

${OBJECTDIR}/src/protonatable-bead.o: src/protonatable-bead.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/protonatable-bead.o src/protonatable-bead.cpp

${OBJECTDIR}/src/protonation-site-catalog.o: src/protonation-site-catalog.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Wall -Iinclude -I../cpputil/include -std=c++14 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/protonation-site-catalog.o src/protonation-site-catalog.cpp

# Subprojects
.build-subprojects:
	cd ../cpputil && ${MAKE}  -f Makefile CONF=Release

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} -r ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcpputil.so
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libparticles.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:
	cd ../cpputil && ${MAKE}  -f Makefile CONF=Release clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
