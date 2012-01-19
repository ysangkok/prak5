////////////////////////////////////////////////////////////////
// GdI 3 - WS 11/12 - Praktikum 2
////////////////////////////////////////////////////////////////
// Bitte tragen Sie hier die Namen und Matrikelnummern Ihrer
// ein (maximal 3).
//
// Bsp.:
// Boris Baumstumpf 0999999133
// Axel Axt 12345678910
// Bruno Schneewittchen 66666336
////////////////////////////////////////////////////////////////

#include <stdexcept>
#include <CL/opencl.h>
#include <string>
#include <fstream>
#include <iostream>

#include <cstdlib>

#include "oclutil.h"


void OCL::setupOcl()
{
	cl_int oclError;

	// get the platform, usually it's only one
	oclError = clGetPlatformIDs(1, &oclPlatform, 0);

	if (oclError != CL_SUCCESS) {
		std::cout << oclErrorString(oclError);
		exit(oclError);
	}

	// query all devices and get the first one
	oclError = clGetDeviceIDs(oclPlatform, gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &oclDevice, 0);
	
	if (oclError != CL_SUCCESS) {
		std::cout << oclErrorString(oclError) << std::endl;
		exit(oclError);
	}

        cl_ulong size;
        clGetDeviceInfo(oclDevice, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
        std::cout << "Local memory size: " << size << std::endl;

	// create a computing context for the device
	oclContext = clCreateContext(0, 1, &oclDevice, 0, 0, &oclError);
	
	if (!oclContext) {
		std::cout << "context not created: " << oclErrorString(oclError) << std::endl;
		exit(oclError);
	}	
	
	// create the command queue for the device and context
	oclCmdQueue = clCreateCommandQueue(oclContext, oclDevice, 0, &oclError);
	
	if (oclError != CL_SUCCESS) {
		std::cout << oclErrorString(oclError) << std::endl;
		exit(oclError);
	}
	
	// load the kernel's source code and create a program
	const char * program = loadProgram().c_str();
	oclProgram = clCreateProgramWithSource(oclContext, 1, &program, 0, &oclError);
	
	if (oclError != CL_SUCCESS) {
		std::cout << oclErrorString(oclError) << std::endl;
		exit(oclError);
	}
	
	// TODO if necessary add optimization options here
	const char * options = "-Werror -g -O0";
	
	// create the kernel with the specified options
	oclError = clBuildProgram(oclProgram, 1, &oclDevice, options, 0, 0);
	
	if (oclError != CL_SUCCESS) {
		std::cout << oclErrorString(oclError) << std::endl;

		{ //if (oclError == CL_BUILD_PROGRAM_FAILURE) {
			size_t length;
			oclError = clGetProgramBuildInfo(oclProgram,
				oclDevice,
				CL_PROGRAM_BUILD_LOG,
				0,
				NULL,
				&length);
			if(oclError != CL_SUCCESS)
				throw std::runtime_error("Can't get program build info (clGetProgramBuildInfo)");
	
			char* buffer = (char*)malloc(length);
			oclError = clGetProgramBuildInfo(oclProgram,
				oclDevice,
				CL_PROGRAM_BUILD_LOG,
				length,
				buffer,
				NULL);
	
			throw std::runtime_error(buffer);
		}

		exit(oclError);
	}
	
	// finally create the kernel
	oclKernel = clCreateKernel(oclProgram, KERNEL_NAME, &oclError);
	
	if (oclError != CL_SUCCESS) {
		std::cout << oclErrorString(oclError) << std::endl;
		exit(oclError);
	}
//	std::cout << "kernel created, ready for computation" << std::endl;
}

void OCL::cleanup()
{
	clReleaseKernel(oclKernel);
	clReleaseProgram(oclProgram);
	clReleaseCommandQueue(oclCmdQueue);
	clReleaseContext(oclContext);
}
