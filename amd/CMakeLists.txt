# CudaCluster main cmake file
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# cluster environment using cuda
PROJECT(Praktikum5)

#------------------------------------------------------------------------------
# options and flags
#------------------------------------------------------------------------------
# Set options
SET(CMAKE_VERBOSE_MAKEFILE off)
SET(CMAKE_BUILD_TYPE "Debug")

SET(BUILD_FOR_64_BITS on CACHE BOOL "Build for 64 bit platforms?")

ADD_DEFINITIONS("-Wall -Wextra")

SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -ffast-math -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")
SET(HEADER ../Fasta.h ../Matrix.h ../MIcomputation.h ../oclutil.h ../SequenceSet.h)
SET(SOURCE ../Fasta.cpp ../pr5.cpp ../SequenceSet.cpp ../MIcomputation.cpp ../oclutil.cpp)

#------------------------------------------------------------------------------
# options and flags
#------------------------------------------------------------------------------
ADD_EXECUTABLE(Praktikum5 ${HEADER} ${SOURCE})

	INCLUDE_DIRECTORIES("/opt/AMDAPP/include")
	find_library(OpenCL OpenCL /usr/lib)
        target_link_libraries(Praktikum5 OpenCL)

