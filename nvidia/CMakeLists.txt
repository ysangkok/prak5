# CudaCluster main cmake file
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# cluster environment using cuda
PROJECT(Praktikum5)

#------------------------------------------------------------------------------
# options and flags
#------------------------------------------------------------------------------
# Set options
SET(CMAKE_VERBOSE_MAKEFILE off)
SET(CMAKE_BUILD_TYPE "Release")

SET(BUILD_FOR_64_BITS on CACHE BOOL "Build for 64 bit platforms?")

SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -march=native -mtune=native")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")
SET(HEADER ../Fasta.h ../Matrix.h ../MIcomputation.h ../oclutil.h ../SequenceSet.h)
SET(SOURCE ../Fasta.cpp ../pr5.cpp ../SequenceSet.cpp ../MIcomputation.cpp ../oclutil.cpp)

# Global include path
INCLUDE_DIRECTORIES("/usr/local/cuda_sdk/OpenCL/common/inc")


#------------------------------------------------------------------------------
# options and flags
#------------------------------------------------------------------------------
ADD_EXECUTABLE(Praktikum5 ${HEADER} ${SOURCE})

TARGET_LINK_LIBRARIES(Praktikum5 "/usr/lib64/libOpenCL.so")
