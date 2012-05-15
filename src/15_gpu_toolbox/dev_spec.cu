/* dev_spec.cu*/

/*
 * Copyright (C) 2008-2012 ABINIT Group (MMancini,FDahm)
 * this file is distributed under the terms of the
 * gnu general public license, see ~abinit/COPYING
 * or http://www.gnu.org/copyleft/gpl.txt.
 *
 */

#include <stdio.h>
#include <string.h>
#include "gpu_four_header.h"

/*=========================================================================*/
/*________________________ GPU_function called by HOST_____________________*/
/*=========================================================================*/
// display CUDA device info
static __host__ void  prt_dev_info()
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  for (int dev = 0; dev < deviceCount; ++dev)
    {
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, dev);
      printf("\n___________________________________________________________________\n");
      printf(  "__________  Graphic Card Properties  ______________________________\n");
      printf("\n  Device %d: \"%s\"\n", dev, deviceProp.name);
      printf("  Revision number:                               %d.%d\n", deviceProp.major,deviceProp.minor);
      printf("  Total amount of global memory:                 %3.1f Mbytes\n", deviceProp.totalGlobalMem/1048576.);
      printf("  Clock rate:                                    %3.1f GHz\n", deviceProp.clockRate/1000000.);
      printf("  Max GFLOP:                                     %d GFP\n", 8*deviceProp.multiProcessorCount * deviceProp.clockRate/1000000);
      printf("  Total amount of constant memory:               %d bytes\n",(int) deviceProp.totalConstMem);
      printf("  Total amount of shared memory per block:       %d bytes\n",(int) deviceProp.sharedMemPerBlock);
      printf("  Total number of registers available per block: %d\n", deviceProp.regsPerBlock);
      printf("___________________________________________________________________\n");
      fflush(stdout);
      if( (int) deviceProp.totalConstMem<0) break;
      //if(deviceProp.major==9999){printf("EXIT: PROBLEM WITH AVAILABLE DEVICES \n");exit(0);}
    }
}


// Explicit Cuda Error ---------------------
__host__  void
check_err(int line ){
/* cuda check errors */
  cudaError_t cudaError;
  cudaError = cudaGetLastError();
  if(cudaError != cudaSuccess)
    { fprintf(stderr, "CUDA Runtime API Error reported : %s %d\n", cudaGetErrorString(cudaError),line);
      exit(EXIT_FAILURE);
    }
  return;
}


// Gives the number of GPU devices ---------
extern "C" __host__
void get_gpu_ndev_(int* ndevice)
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  *ndevice = deviceCount;

  return;
}

// Gives the max memory available for a GPU device ---------
extern "C" __host__
void get_gpu_max_mem_(int* device, float* max_mem)
{
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, *device);
  *max_mem = (float) deviceProp.totalGlobalMem;
   return;
}


// Set the device if it exists   -----------------
extern "C" __host__
void set_dev_(int* gpudevice)
{
 if(*gpudevice >-1){
   cudaError_t cudaError;
   int deviceCount;
   cudaGetDeviceCount(&deviceCount);
   if(deviceCount>*gpudevice){
     cudaSetDevice(*gpudevice);
     cudaError = cudaGetLastError();
     if(cudaError != cudaSuccess){
       fprintf(stderr, "CUDA Runtime API Error reported : %s\n", cudaGetErrorString(cudaError));
       fflush(stderr);
       exit(1);
     }
   }
   else *gpudevice=-1;
 }
  return;
}


// Unset the devices  -----------------
extern "C"  __host__
void unset_dev_()
{
#if defined HAVE_GPU_CUDA3
  cudaThreadExit();
#else
  cudaDeviceReset();
#endif
  return;
}


// Get context  -----------------------
extern "C"  __host__
void check_context_(int *res,char *message)
{
  *res=1;
  cudaError_t state=cudaFree(0);
  if (state!=cudaSuccess){
    sprintf(message,"Unable to initialize a Cuda context: %s \n",cudaGetErrorString(state));
    *res=0;
    unset_dev_();
  }
}


// Get info from device  --------------
extern "C" __host__
void  get_dev_info_(int* device,
		    char* name,
		    int* lenname,
		    int vers[2],
		    float* globalmem,
		    float* clockrate,
		    int* gflop,
		    int* constmem,
		    int* sharemem,
		    int* regist
		    )
{
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, *device);
  strcpy(name,deviceProp.name);
  *lenname = strlen( name );
  vers[0] = deviceProp.major;
  vers[1] = deviceProp.minor;
  *globalmem = deviceProp.totalGlobalMem/1048576.;
  *clockrate = deviceProp.clockRate/1000000.;
  *gflop = 8*deviceProp.multiProcessorCount * deviceProp.clockRate/1000000;
  *constmem = deviceProp.totalConstMem;
  *sharemem =  deviceProp.sharedMemPerBlock;
  *regist = deviceProp.regsPerBlock;
}


extern "C"  __host__
void c_get_ndevice_(int* ndev)
{
  *ndev=0;
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  for (int idev = 0; idev < deviceCount; ++idev)
    {
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, idev);
      //We check that no device is in "emu" mode
      if( deviceProp.major != 9999 ) {
#if defined HAVE_GPU_CUDA_DP
      //We check that double precision is available, c.c. >= 1.3 )
      if( (deviceProp.major>1)||(deviceProp.minor>2) )
#endif
	*ndev+=1;
      }
    }
}


/***************************************************************/
/*******                                                ********/
/*******      GPU MEMORY MANAGEMENT ROUTINES            ********/
/*******                                                ********/
/***************************************************************/

/*============================================================================*/
/* Print memory information (total amount and free available)                 */
/*============================================================================*/

extern "C" void check_gpu_mem_(){
  size_t free,total;
  cudaMemGetInfo(&free,&total);
  printf("*** GPU memory : Free =>  %4.2fMo   | Total =>  %4.2fMo ***\n",free*1e-6,total*1e-6);
  fflush(stdout);
}

/*============================================================================*/
/* Allocate size byte in gpu memory and returns in gpu_ptr this location      */
/* INPUTS size= size in byte to allocate                                      */
/* OUTPUT gpu_ptr= C_PTR on gpu memory location that has been allocated       */
/*============================================================================*/

extern "C" void alloc_on_gpu_(void **gpu_ptr,int* size){

  if(cudaMalloc(gpu_ptr,*size)!=cudaSuccess){
    printf("ERROR: alloc_on_gpu failed:%s\n",cudaGetErrorString(cudaGetLastError()));
    fflush(stdout);
    leave_new_("COLL");
  }
}

/*============================================================================*/
/* Free memory location pointed by gpu_ptr                                    */
/* OUTPUT gpu_ptr= C_PTR on gpu memory location that has been allocated       */
/* WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled   */
/*            the correct one is in xx_gpu_toolbox/dev_spec.cu                */
/*============================================================================*/

extern "C" void dealloc_on_gpu_(void **gpu_ptr){
  if(*gpu_ptr==NULL)
    return;
  if(cudaFree(*gpu_ptr)!=cudaSuccess){
    printf("ERROR: dealloc_on_gpu failed :%s\n",cudaGetErrorString(cudaGetLastError()));
    fflush(stdout);
    leave_new_("COLL");
  }
  *gpu_ptr=NULL;
}

/*============================================================================*/
/* Copy size byte from  dtab to gpu memory pointed by gpu_ptr                 */
/* INPUTS                                                                     */
/*  size= size in byte to allocate                                            */
/*  dtab = fortran tab to copy                                                */
/* OUTPUT                                                                     */
/*  gpu_ptr= C_PTR on gpu memory location                                     */
/* WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled   */
/*            the correct one is in xx_gpu_toolbox/dev_spec.cu                */
/*============================================================================*/

extern "C" void copy_on_gpu_(void *ptr, void **gpu_ptr,int* size){
  if(cudaMemcpy(*gpu_ptr,ptr,*size,cudaMemcpyHostToDevice)!=cudaSuccess){
    printf("ERROR: copy_on_gpu failed : %s\n",cudaGetErrorString(cudaGetLastError()));
    fflush(stdout);
    leave_new_("COLL");
  }
}

/*============================================================================*/
/* Copy size byte from gpu memory pointed by gpu_ptr to dtab                  */
/* INPUTS                                                                     */
/*  size= size in byte to allocate                                            */
/*  gpu_ptr= C_PTR on gpu memory location that has been allocated             */
/* OUTPUT                                                                     */
/*  dtab = fortran tab which will contains data                               */
/*============================================================================*/

extern "C" void copy_from_gpu_(void *ptr,void **gpu_ptr,int* size){
  if(cudaMemcpy(ptr,*gpu_ptr,*size,cudaMemcpyDeviceToHost)!=cudaSuccess){
    printf("ERROR: copy_from_gpu failed : %s\n",cudaGetErrorString(cudaGetLastError()));
    fflush(stdout);
    leave_new_("COLL");
  }
}
